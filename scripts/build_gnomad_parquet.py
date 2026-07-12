#!/usr/bin/env python3
"""Build a gnomAD frequency Parquet from scratch, straight from the internet.

Reproduces the genomic-lakehouse gnomAD store: for each chromosome it STREAMS the
gnomAD sites VCF (never storing the raw ~150-200 GB — a sequential scan is bandwidth-
bound, not the per-variant-latency trap), extracts only the fields the ACMG engine
cites, and writes a Hive-partitioned Parquet (``chrom=chrN/data.parquet``). DuckDB
then answers a whole exome's frequencies in one vectorised join in ~seconds, offline.

Schema (per row): chrom, pos, ref, alt, filter, af, af_grpmax, ac, an, nhomalt,
faf95, grpmax_pop, and per-ancestry AFs af_afr/af_amr/af_asj/af_eas/af_fin/af_mid/
af_nfe/af_sas/af_ami/af_remaining. ``faf95`` is gnomAD's grpmax filtering AF (ClinGen
BS1/BA1) — the joint release carries it (``fafmax_faf95_max_joint``), an improvement
over an af_grpmax-only table; the per-population AFs reach parity with a lakehouse store.

    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_parquet.py --chroms 21
    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_parquet.py           # all
    python3 scripts/build_gnomad_parquet.py --src /Volumes/DISK/gnomad_joint     # local VCFs

Needs ``bcftools`` (stream/extract) and ``duckdb`` (write Parquet). ~1 GB output for
the full v4.1 joint (29.6M variants, incl. per-population AFs); the raw VCFs are never kept.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

# preset -> (per-chrom URL template, ordered [(out_col, INFO field)]).
# Fixed leading columns are always CHROM,POS,REF,ALT,FILTER.
_GCS = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
# gnomAD v4 genetic-ancestry groups -> per-population AF columns (af_afr, af_amr, ...),
# for ancestry-aware interpretation (a variant common in one group but not the grpmax).
# 'raw' is the unfiltered AF, not an ancestry group, so it is excluded. The exomes-only
# release has no 'ami' (Amish) group; the joint release carries it (via genomes).
_POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]
_POPS_EXOMES = [p for p in _POPS if p != "ami"]


def _pops(info_prefix: str, pops):
    """[(af_afr, <prefix>_afr), ...] — joint fields are AF_joint_<pop>, exomes AF_<pop>."""
    return [(f"af_{p}", f"{info_prefix}_{p}") for p in pops]


_PRESETS = {
    "joint": (
        f"{_GCS}/joint/gnomad.joint.v4.1.sites.{{chrom}}.vcf.bgz",
        [("af", "AF_joint"), ("af_grpmax", "AF_grpmax_joint"), ("ac", "AC_joint"),
         ("an", "AN_joint"), ("nhomalt", "nhomalt_joint"),
         ("faf95", "fafmax_faf95_max_joint"), ("grpmax_pop", "grpmax_joint")]
        + _pops("AF_joint", _POPS),
    ),
    "exomes": (
        f"{_GCS}/exomes/gnomad.exomes.v4.1.sites.{{chrom}}.vcf.bgz",
        [("af", "AF"), ("af_grpmax", "AF_grpmax"), ("ac", "AC"), ("an", "AN"),
         ("nhomalt", "nhomalt"), ("faf95", "fafmax_faf95_max"), ("grpmax_pop", "grpmax")]
        + _pops("AF", _POPS_EXOMES),
    ),
}
_LEAD = ["chrom", "pos", "ref", "alt", "filter"]
_NUMERIC = {"pos": "INTEGER", "af": "DOUBLE", "af_grpmax": "DOUBLE", "ac": "BIGINT",
            "an": "BIGINT", "nhomalt": "BIGINT", "faf95": "DOUBLE",
            **{f"af_{p}": "DOUBLE" for p in _POPS}}


def _which(tool: str) -> bool:
    from shutil import which
    return which(tool) is not None


def _chroms(spec: str) -> list[str]:
    """'21' / '1-22,X,Y' / '1,2,3' -> ['21'] / ['1'..'22','X','Y']."""
    out: list[str] = []
    for part in spec.split(","):
        part = part.strip()
        if "-" in part and part.replace("-", "").isdigit():
            a, b = part.split("-")
            out += [str(i) for i in range(int(a), int(b) + 1)]
        elif part:
            out.append(part)
    return out


def _source(url_tmpl: str, chrom: str, src: str | None) -> str:
    label = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    if src:
        return str(Path(src) / Path(url_tmpl.format(chrom=label)).name)
    return url_tmpl.format(chrom=label)


def _bcftools_fmt(fields) -> str:
    cols = ["%CHROM", "%POS", "%REF", "%ALT", "%FILTER"]
    cols += [f"%INFO/{info}" for _out, info in fields]
    return "\\t".join(cols) + "\\n"


def _copy_sql(tsv: Path, out_parquet: Path, cols: list[str]) -> str:
    read_cols = ", ".join(f"'{c}': 'VARCHAR'" for c in cols)
    sel = []
    for c in cols:
        if c in _NUMERIC:
            sel.append(f"TRY_CAST(NULLIF({c}, '.') AS {_NUMERIC[c]}) AS {c}")
        else:
            sel.append(f"NULLIF({c}, '.') AS {c}")
    return (
        f"COPY (SELECT {', '.join(sel)} FROM read_csv('{tsv}', delim='\\t', "
        f"header=false, columns={{{read_cols}}})) "
        f"TO '{out_parquet}' (FORMAT PARQUET, COMPRESSION ZSTD);"
    )


def _download(url_tmpl: str, chrom: str, dl_dir: str) -> tuple[str, list[Path]]:
    """Resumably download a chrom's bgz + .tbi to dl_dir (curl -C - --retry). Returns
    (local bgz path, [files to delete]). Local processing then never stalls on a remote
    read (the streaming failure mode) and needs only ~one chrom of scratch at a time."""
    label = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    url = url_tmpl.format(chrom=label)
    Path(dl_dir).mkdir(parents=True, exist_ok=True)
    bgz = Path(dl_dir) / Path(url).name
    tbi = Path(str(bgz) + ".tbi")
    for u, dst in [(url + ".tbi", tbi), (url, bgz)]:   # index first (small), then data
        print(f"[{label}] downloading {dst.name} ...", file=sys.stderr)
        subprocess.run(["curl", "-fsSL", "-C", "-", "--retry", "10", "--retry-delay", "5",
                        "-o", str(dst), u], check=True)
    return str(bgz), [bgz, tbi]


def build_chrom(chrom: str, url_tmpl: str, fields, src: str | None,
                out_dir: Path, region: str | None, bed: str | None = None,
                download_dir: str | None = None) -> int:
    import duckdb
    label = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    cols = _LEAD + [out for out, _info in fields]
    part_dir = out_dir / f"chrom={label}"
    part_dir.mkdir(parents=True, exist_ok=True)
    tsv = part_dir / "_extract.tsv"
    parquet = part_dir / "data.parquet"

    if parquet.exists() and parquet.stat().st_size > 0:   # resume: this chrom is done
        con = duckdb.connect()
        n = con.execute(f"SELECT count(*) FROM read_parquet('{parquet}')").fetchone()[0]
        con.close()
        print(f"[{label}] {n:,} variants (already built — skipped)", file=sys.stderr)
        return n

    cleanup: list[Path] = []
    if download_dir and not src:
        source, cleanup = _download(url_tmpl, chrom, download_dir)
    else:
        source = _source(url_tmpl, chrom, src)

    cmd = ["bcftools", "query", "-f", _bcftools_fmt(fields).replace("\\t", "\t").replace("\\n", "\n")]
    if region:
        cmd += ["-r", region]
    if bed:
        # Restrict to this chrom's target regions (the genomic-lakehouse method: fetch only
        # the exome-panel intervals via the tabix index, not the whole ~17-67 GB file).
        chrom_bed = part_dir / "_regions.bed"
        want = (label, chrom)
        rows = [ln for ln in Path(bed).read_text().splitlines()
                if ln and not ln.startswith(("#", "track", "browser"))
                and ln.split("\t", 1)[0] in want]
        chrom_bed.write_text("\n".join(rows) + "\n")
        if not rows:
            print(f"[{label}] no BED regions for this contig — skipped", file=sys.stderr)
            for f in cleanup:
                f.unlink(missing_ok=True)
            return 0
        cmd += ["-R", str(chrom_bed)]
    cmd.append(source)
    print(f"[{label}] {'extracting from ' + Path(source).name if cleanup or src else 'streaming ' + source} ...",
          file=sys.stderr)
    with open(tsv, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if p.returncode != 0:
        sys.stderr.write(p.stderr.decode("utf-8", "replace")[-500:] + "\n")
        tsv.unlink(missing_ok=True)
        for f in cleanup:                 # keep the (resumable) download for a retry? no —
            f.unlink(missing_ok=True)     # a failed extract is likely a bad/partial file
        return 0

    con = duckdb.connect()
    con.execute(_copy_sql(tsv, parquet, cols))
    n = con.execute(f"SELECT count(*) FROM read_parquet('{parquet}')").fetchone()[0]
    con.close()
    tsv.unlink(missing_ok=True)           # never keep the raw extract
    for f in cleanup:                     # free the ~one-chrom scratch before the next
        f.unlink(missing_ok=True)
    print(f"[{label}] {n:,} variants -> {parquet}", file=sys.stderr)
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--preset", choices=list(_PRESETS), default="joint",
                    help="gnomAD dataset (default: joint = exomes+genomes v4.1).")
    ap.add_argument("--chroms", default="1-22,X,Y", help="e.g. 21 | 1-22,X,Y")
    ap.add_argument("--src", help="dir of local per-chrom gnomAD VCFs (no network).")
    ap.add_argument("--download-dir", help="download each chrom's VCF here (resumable), "
                    "process it locally, then delete it — robust vs streaming (which stalls "
                    "on remote reads) and needs only ~one chrom of scratch at a time. Point "
                    "at a disk with room (e.g. an external SSD).")
    ap.add_argument("--region", help="a single region (e.g. chr21:31659622-31668931) "
                                     "for a quick end-to-end test on one chromosome.")
    ap.add_argument("--bed", help="restrict to an exome-panel BED (the genomic-lakehouse "
                                  "method: fetch only the panel intervals via the tabix "
                                  "index — ~600 MB total, feasible from a laptop). Best "
                                  "with --preset exomes. Produces a panel-scoped store.")
    ap.add_argument("--out", default=str(Path(config.GNOMAD_LOCAL_TABIX).parent / "gnomad_parquet"),
                    help="output partitioned-parquet dir (default: %(default)s).")
    args = ap.parse_args(argv)

    if not args.src and config.offline():
        print("ERROR: needs network to stream gnomAD (or pass --src DIR of local VCFs).\n"
              "  Set VCF2REPORT_ALLOW_NETWORK=1.", file=sys.stderr)
        return 2
    if not _which("bcftools"):
        print("ERROR: bcftools not found (conda install -c bioconda bcftools).", file=sys.stderr)
        return 2
    try:
        import duckdb  # noqa: F401
    except Exception:
        print("ERROR: duckdb not found (pip install duckdb).", file=sys.stderr)
        return 2

    url_tmpl, fields = _PRESETS[args.preset]
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    chroms = _chroms(args.region.split(":")[0] if args.region else args.chroms)

    full_set = {str(i) for i in range(1, 23)} | {"X", "Y"}
    built: dict[str, int] = {}
    total = 0
    for c in chroms:
        n = build_chrom(c, url_tmpl, fields, args.src, out_dir, args.region, args.bed,
                        args.download_dir)
        built[c] = n
        total += n

    ok = {c for c, n in built.items() if n > 0}
    # Declare 'full' ONLY for a whole-genome build where every requested chromosome
    # streamed successfully and there was no --region/--bed slice — so the client may
    # assert a variant absent from these contigs. Anything else is 'partial': the client
    # then never fabricates an absence off it (a region test, a panel, or a truncated
    # build). A --bed store is panel-scoped -> always partial (absence off-panel is unknown).
    is_full = (not args.region) and (not args.bed) and (set(chroms) >= full_set) and (ok >= set(chroms))
    _pref = lambda c: c if str(c).lower().startswith("chr") else f"chr{c}"
    meta = {"mode": "full" if is_full else "partial",
            "contigs": sorted(_pref(c) for c in ok),
            "preset": args.preset, "rows": total,
            "source": f"gnomAD v4.1 {args.preset}"}
    (out_dir / "_meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    print(f"\nDone ({meta['mode']}): {total:,} variants across {len(ok)} chrom(s) -> {out_dir}",
          file=sys.stderr)
    if not is_full and not args.region and not args.bed:
        _warn("mode=partial: some chromosomes did not build — the client will NOT assert "
              "absence off this store (falls back instead). Re-run to complete for full mode.")
    print(f"Point vcf2report at it: VCF2REPORT_GNOMAD_PARQUET={out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
