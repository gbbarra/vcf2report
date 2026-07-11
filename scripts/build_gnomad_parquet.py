#!/usr/bin/env python3
"""Build a gnomAD frequency Parquet from scratch, straight from the internet.

Reproduces the genomic-lakehouse gnomAD store: for each chromosome it STREAMS the
gnomAD sites VCF (never storing the raw ~150-200 GB — a sequential scan is bandwidth-
bound, not the per-variant-latency trap), extracts only the fields the ACMG engine
cites, and writes a Hive-partitioned Parquet (``chrom=chrN/data.parquet``). DuckDB
then answers a whole exome's frequencies in one vectorised join in ~seconds, offline.

Schema (per row): chrom, pos, ref, alt, filter, af, af_grpmax, ac, an, nhomalt,
faf95, grpmax_pop. ``faf95`` is gnomAD's grpmax filtering AF (ClinGen BS1/BA1) — the
joint release carries it (``fafmax_faf95_max_joint``), an improvement over an
af_grpmax-only table.

    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_parquet.py --chroms 21
    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_parquet.py           # all
    python3 scripts/build_gnomad_parquet.py --src /Volumes/DISK/gnomad_joint     # local VCFs

Needs ``bcftools`` (stream/extract) and ``duckdb`` (write Parquet). ~786 MB output for
the full v4.1 joint (29.6M variants); the raw VCFs are never kept.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

# preset -> (per-chrom URL template, ordered [(out_col, INFO field)]).
# Fixed leading columns are always CHROM,POS,REF,ALT,FILTER.
_GCS = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
_PRESETS = {
    "joint": (
        f"{_GCS}/joint/gnomad.joint.v4.1.sites.{{chrom}}.vcf.bgz",
        [("af", "AF_joint"), ("af_grpmax", "AF_grpmax_joint"), ("ac", "AC_joint"),
         ("an", "AN_joint"), ("nhomalt", "nhomalt_joint"),
         ("faf95", "fafmax_faf95_max_joint"), ("grpmax_pop", "grpmax_joint")],
    ),
    "exomes": (
        f"{_GCS}/exomes/gnomad.exomes.v4.1.sites.{{chrom}}.vcf.bgz",
        [("af", "AF"), ("af_grpmax", "AF_grpmax"), ("ac", "AC"), ("an", "AN"),
         ("nhomalt", "nhomalt"), ("faf95", "fafmax_faf95_max"), ("grpmax_pop", "grpmax")],
    ),
}
_LEAD = ["chrom", "pos", "ref", "alt", "filter"]
_NUMERIC = {"pos": "INTEGER", "af": "DOUBLE", "af_grpmax": "DOUBLE", "ac": "BIGINT",
            "an": "BIGINT", "nhomalt": "BIGINT", "faf95": "DOUBLE"}


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


def build_chrom(chrom: str, url_tmpl: str, fields, src: str | None,
                out_dir: Path, region: str | None) -> int:
    import duckdb
    label = chrom if chrom.lower().startswith("chr") else f"chr{chrom}"
    source = _source(url_tmpl, chrom, src)
    cols = _LEAD + [out for out, _info in fields]
    part_dir = out_dir / f"chrom={label}"
    part_dir.mkdir(parents=True, exist_ok=True)
    tsv = part_dir / "_extract.tsv"
    parquet = part_dir / "data.parquet"

    cmd = ["bcftools", "query", "-f", _bcftools_fmt(fields).replace("\\t", "\t").replace("\\n", "\n")]
    if region:
        cmd += ["-r", region]
    cmd.append(source)
    print(f"[{label}] streaming {source} ...", file=sys.stderr)
    with open(tsv, "w") as fh:
        p = subprocess.run(cmd, stdout=fh, stderr=subprocess.PIPE)
    if p.returncode != 0:
        sys.stderr.write(p.stderr.decode("utf-8", "replace")[-500:] + "\n")
        tsv.unlink(missing_ok=True)
        return 0

    con = duckdb.connect()
    con.execute(_copy_sql(tsv, parquet, cols))
    n = con.execute(f"SELECT count(*) FROM read_parquet('{parquet}')").fetchone()[0]
    con.close()
    tsv.unlink(missing_ok=True)   # never keep the raw extract
    print(f"[{label}] {n:,} variants -> {parquet}", file=sys.stderr)
    return n


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--preset", choices=list(_PRESETS), default="joint",
                    help="gnomAD dataset (default: joint = exomes+genomes v4.1).")
    ap.add_argument("--chroms", default="1-22,X,Y", help="e.g. 21 | 1-22,X,Y")
    ap.add_argument("--src", help="dir of local per-chrom gnomAD VCFs (no network).")
    ap.add_argument("--region", help="a single region (e.g. chr21:31659622-31668931) "
                                     "for a quick end-to-end test on one chromosome.")
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

    total = 0
    for c in chroms:
        total += build_chrom(c, url_tmpl, fields, args.src, out_dir, args.region)

    (out_dir / "_SUCCESS").write_text(f"preset={args.preset} chroms={','.join(chroms)} rows={total}\n")
    print(f"\nDone: {total:,} variants across {len(chroms)} chrom(s) -> {out_dir}", file=sys.stderr)
    print(f"Point vcf2report at it: VCF2REPORT_GNOMAD_PARQUET={out_dir}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
