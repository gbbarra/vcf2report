#!/usr/bin/env python3
"""Build the reduced, offline, tabix-indexed gnomAD frequency table.

Produces the bgzipped TSV (+ .tbi) that ``src/vcf2report/annotate/gnomad_local.py``
reads. Every reduced record is computed by REUSING ``gnomad_remote``'s per-record
reduction (grpmax AF/AC/AN/hom + ClinGen filtering AF ``faf95``, higher grpmax of the
exomes vs genomes callsets), so a local hit equals what the live remote path would
return — the local table just makes it offline and instant.

Schema (locked — the client depends on it EXACTLY; tab-separated, sorted by chrom then
pos, ``chr`` stripped, ``#`` header comment first):

    #chrom  pos  ref  alt  af  ac  an  hom  faf95  pop

Modes
-----
  --from-vcf VCF   (default: the bundled demo VCF) one row per input variant via
                   ``gnomad_remote.query`` — covers exactly this VCF's sites. Partial.
  --bed BED        region-tabix gnomAD for each BED interval; one reduced row per
                   gnomAD variant in the region. Partial.
  --full [--src D] every variant of chr1..22,X,Y. With --src, reads local per-chrom
                   sites files; otherwise streams the public GCS bucket (~150-200 GB).
                   Full — the only build allowed to assert a true absence on a miss.

Examples
--------
    # tiny per-VCF table (the common case) — needs network for the gnomAD lookups
    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_gnomad_local.py --from-vcf my.vcf

    # a panel of regions
    VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_gnomad_local.py --bed panel.bed

    # full genome-wide table from local sites files (no network)
    python scripts/build_gnomad_local.py --full --src /Volumes/DATA/gnomad_v4.1
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import config                        # noqa: E402
from vcf2report.annotate import gnomad_remote        # noqa: E402
from vcf2report.vcf.parse import parse_vcf           # noqa: E402
from vcf2report.models import Variant                # noqa: E402

HEADER = "#chrom\tpos\tref\talt\taf\tac\tan\thom\tfaf95\tpop"
RELEASE = gnomad_remote.RELEASE
SOURCE = f"gnomAD v{RELEASE} exomes+genomes"
CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y"]


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _warn(msg: str) -> None:
    print(f"WARNING: {msg}", file=sys.stderr)


def _norm_chrom(c: str) -> str:
    c = str(c)
    return c[3:] if c.lower().startswith("chr") else c


def _prefixed(c: str) -> str:
    """Contig name as it appears in gnomAD files (always ``chr``-prefixed)."""
    c = str(c)
    return c if c.lower().startswith("chr") else f"chr{c}"


def _fmt(v) -> str:
    """Render a reduced field for the TSV: None -> '', float -> round-tripping repr."""
    if v is None:
        return ""
    if isinstance(v, float):
        return repr(v)
    return str(v)


def _row(chrom, pos, ref, alt, r: dict) -> str:
    return "\t".join((
        _norm_chrom(chrom),
        str(int(pos)),
        str(ref),
        str(alt),
        _fmt(r.get("af")),
        _fmt(r.get("ac")),
        _fmt(r.get("an")),
        _fmt(r.get("hom")),
        _fmt(r.get("faf95")),
        (r.get("pop") or ""),
    ))


def _to_float(s: str):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.1f}{unit}"
        n /= 1024
    return f"{n:.1f}PB"


# ---------------------------------------------------------------------------
# Source opening / record reduction (reuse gnomad_remote so local == remote)
# ---------------------------------------------------------------------------
def _open_source(kind: str, chrom_pref: str, src: str | None):
    """Return a VariantFile for (kind, chrom): local sites file under ``src`` if
    given, else the cached remote handle from ``gnomad_remote``. None if missing."""
    if src:
        import pysam
        p = Path(src) / f"gnomad.{kind}.v{RELEASE}.sites.{chrom_pref}.vcf.bgz"
        if not p.exists():
            _warn(f"missing local sites file: {p}")
            return None
        try:
            return pysam.VariantFile(str(p))
        except Exception as e:                       # unreadable / bad index
            _warn(f"cannot open {p}: {e}")
            return None
    return gnomad_remote._open(kind, chrom_pref)     # cached, remote GCS URL


def _emit_records(chrom_norm, chrom_pref, start, end, src, raw_fh) -> int:
    """Reduce every gnomAD variant of both callsets in a chrom (start/end None) or a
    region into raw rows, writing them UNMERGED — the post-sort collapse keeps the
    higher-grpmax duplicate, matching gnomad_remote.query's exomes-vs-genomes merge."""
    n = 0
    for kind in ("exomes", "genomes"):
        vf = _open_source(kind, chrom_pref, src)
        if vf is None:
            continue
        try:
            recs = vf.fetch(chrom_pref) if start is None else vf.fetch(chrom_pref, start, end)
        except Exception as e:
            _warn(f"fetch failed for {kind} {chrom_pref}"
                  f"{'' if start is None else f':{start}-{end}'}: {e}")
            continue
        try:
            for rec in recs:
                alts = rec.alts or ()
                if not alts:
                    continue
                cand = gnomad_remote._best_from_record(rec)  # grpmax AF/faf95 reduction
                if cand is None:
                    continue
                raw_fh.write(_row(chrom_norm, rec.pos, rec.ref, alts[0], cand) + "\n")
                n += 1
        except Exception as e:                       # flaky remote read mid-stream
            _warn(f"stream interrupted for {kind} {chrom_pref}: {e}")
            continue
    return n


# ---------------------------------------------------------------------------
# Generators (one per mode). Each streams raw rows into ``raw_fh``.
# ---------------------------------------------------------------------------
def _query_with_retry(v: Variant, attempts: int = 3):
    """gnomad_remote.query with retry: None is a transport failure (retry+backoff),
    the absent dict (af 0.0) is a real covered-but-absent answer (returned)."""
    for a in range(1, attempts + 1):
        r = gnomad_remote.query(v)
        if r is not None:
            return r
        if a < attempts:
            time.sleep(0.5 * a)
    return None


def _gen_from_vcf(vcf_path: Path, raw_fh, progress_every: int) -> tuple[int, int]:
    variants, build, _hdr = parse_vcf(vcf_path)
    if build and build != config.GENOME_BUILD:
        _warn(f"VCF build {build!r} != {config.GENOME_BUILD}; gnomAD v{RELEASE} is "
              f"{config.GENOME_BUILD} — coordinates may not match.")
    seen: set[str] = set()
    uniq = [v for v in variants if not (v.key in seen or seen.add(v.key))]
    print(f"Querying gnomAD for {len(uniq)} unique sites "
          f"({len(variants)} variants in {vcf_path})...", file=sys.stderr)
    written = skipped = 0
    for i, v in enumerate(uniq, 1):
        r = _query_with_retry(v)
        if r is None:                                # transport failure -> skip, no fabrication
            skipped += 1
            _warn(f"skipping {v.key}: gnomAD lookup failed after retries")
            continue
        raw_fh.write(_row(v.chrom, v.pos, v.ref, v.alt, r) + "\n")
        written += 1
        if i % progress_every == 0:
            print(f"  ...{i}/{len(uniq)} sites ({written} rows, {skipped} skipped)",
                  file=sys.stderr)
    return written, skipped


def _parse_bed(bed_path: Path) -> list[tuple[str, int, int]]:
    regions: list[tuple[str, int, int]] = []
    for ln in bed_path.read_text().splitlines():
        ln = ln.strip()
        if not ln or ln.startswith(("#", "track", "browser")):
            continue
        f = ln.split("\t") if "\t" in ln else ln.split()
        if len(f) < 3:
            continue
        try:
            regions.append((f[0], int(f[1]), int(f[2])))   # BED: 0-based, half-open
        except ValueError:
            _warn(f"skipping malformed BED line: {ln}")
    return regions


def _gen_bed(bed_path: Path, src: str | None, raw_fh, progress_every: int) -> tuple[int, int]:
    regions = _parse_bed(bed_path)
    print(f"Reducing gnomAD over {len(regions)} BED region(s)...", file=sys.stderr)
    written = 0
    for j, (chrom, start, end) in enumerate(regions, 1):
        written += _emit_records(_norm_chrom(chrom), _prefixed(chrom), start, end, src, raw_fh)
        if j % progress_every == 0:
            print(f"  ...{j}/{len(regions)} regions ({written} rows)", file=sys.stderr)
    return written, 0


def _gen_full(src: str | None, raw_fh, progress_every: int) -> tuple[int, int]:
    written = 0
    for chrom in CHROMS:
        before = written
        written += _emit_records(chrom, _prefixed(chrom), None, None, src, raw_fh)
        print(f"  chr{chrom}: {written - before} variants "
              f"({written} total)", file=sys.stderr)
    return written, 0


# ---------------------------------------------------------------------------
# Sort + collapse + index
# ---------------------------------------------------------------------------
def _sort_file(raw_path: Path, sorted_path: Path) -> None:
    """Sort raw rows by (chrom, pos-numeric, ref, alt). Tries the streaming system
    ``sort`` (handles the huge --full case); falls back to an in-memory Python sort."""
    env = dict(os.environ, LC_ALL="C")
    try:
        subprocess.run(
            ["sort", "-t", "\t", "-k1,1", "-k2,2n", "-k3,3", "-k4,4",
             "-o", str(sorted_path), str(raw_path)],
            check=True, env=env)
        return
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        _warn(f"system sort unavailable ({e}); sorting in Python (needs RAM).")
    rows = [ln for ln in raw_path.read_text().splitlines() if ln]

    def _key(ln: str):
        f = ln.split("\t")
        try:
            return (f[0], int(f[1]), f[2], f[3])
        except (IndexError, ValueError):
            return (f[0] if f else "", 0, "", "")

    rows.sort(key=_key)
    sorted_path.write_text(("\n".join(rows) + "\n") if rows else "")


def _collapse(sorted_path: Path, final_path: Path) -> int:
    """Write HEADER + the sorted rows, collapsing consecutive duplicate
    (chrom,pos,ref,alt) keys to the single higher-``af`` row (exomes vs genomes)."""
    kept = 0
    cur_key = None
    best_line = None
    best_af = float("-inf")
    with open(sorted_path) as inp, open(final_path, "w") as out:
        out.write(HEADER + "\n")
        for line in inp:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) < 4:
                continue
            key = (f[0], f[1], f[2], f[3])
            af = _to_float(f[4]) if len(f) > 4 else None
            afv = af if af is not None else float("-inf")
            if key != cur_key:
                if best_line is not None:
                    out.write(best_line + "\n")
                    kept += 1
                cur_key, best_line, best_af = key, line, afv
            elif afv > best_af:
                best_line, best_af = line, afv
        if best_line is not None:
            out.write(best_line + "\n")
            kept += 1
    return kept


def _index(final_path: Path, out: Path) -> None:
    import pysam
    pysam.tabix_compress(str(final_path), str(out), force=True)
    pysam.tabix_index(str(out), seq_col=0, start_col=1, end_col=1,
                      meta_char="#", line_skip=0, force=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build the reduced offline gnomAD frequency table "
                    "(bgzipped TSV + .tbi) read by annotate/gnomad_local.py.")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--from-vcf", metavar="VCF",
                   help="build a per-VCF table (default: the bundled demo VCF).")
    g.add_argument("--bed", metavar="BED",
                   help="build a table covering the given BED regions.")
    g.add_argument("--full", action="store_true",
                   help="build a genome-wide table (chr1..22,X,Y). Heavy.")
    p.add_argument("--src", metavar="DIR",
                   help="[--full only] directory of local gnomAD per-chrom sites "
                        "files (gnomad.<kind>.v%s.sites.chrN.vcf.bgz); no network." % RELEASE)
    p.add_argument("--out", metavar="PATH", default=str(config.GNOMAD_LOCAL_TABIX),
                   help="output .tsv.gz path (default: %(default)s).")
    p.add_argument("--progress-every", type=int, default=None,
                   help="progress interval (default: adaptive per mode).")
    return p


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)

    # Resolve mode + input.
    if args.full:
        mode, built_from = "full", (f"--full --src {args.src}" if args.src else "--full (remote)")
    elif args.bed:
        mode, built_from = "bed", f"--bed {args.bed}"
    else:
        vcf = Path(args.from_vcf) if args.from_vcf else config.SAMPLE_VCF
        mode, built_from = "from-vcf", f"--from-vcf {vcf}"
    if args.src and not args.full:
        print("ERROR: --src is only valid with --full.", file=sys.stderr)
        return 2

    table_mode = "full" if mode == "full" else "partial"   # sidecar "mode" field

    # Network gating: everything but --full --src reaches gnomAD over the network.
    needs_net = not (mode == "full" and args.src)
    if needs_net and config.offline():
        print("ERROR: this build needs network access to reach gnomAD, but egress is "
              "disabled.\n  Set VCF2REPORT_ALLOW_NETWORK=1 (and ensure OFFLINE is unset).",
              file=sys.stderr)
        return 2

    # Verify pysam is importable up front (needed for tabix compress/index).
    try:
        import pysam  # noqa: F401
    except Exception:
        print("ERROR: pysam is required (tabix compress/index). pip install pysam",
              file=sys.stderr)
        return 2

    if mode == "full":
        print("=" * 72, file=sys.stderr)
        print("WARNING: --full is a HEAVY build.", file=sys.stderr)
        if args.src:
            print(f"  Reading local gnomAD v{RELEASE} sites files from: {args.src}",
                  file=sys.stderr)
        else:
            print(f"  Streaming gnomAD v{RELEASE} exomes+genomes from GCS: this reads on"
                  " the order of ~150-200 GB over the network.", file=sys.stderr)
        print("  The output table is genome-wide; ensure the --out volume has room.",
              file=sys.stderr)
        print("=" * 72, file=sys.stderr)

    progress_every = args.progress_every or {"from-vcf": 25, "bed": 500, "full": 1}[mode]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    raw = out.parent / (out.name + ".build.raw.tmp")
    srt = out.parent / (out.name + ".build.sorted.tmp")
    fin = out.parent / (out.name + ".build.final.tmp")

    written = skipped = 0
    try:
        with open(raw, "w") as raw_fh:
            if mode == "from-vcf":
                if not vcf.exists():
                    print(f"ERROR: VCF not found: {vcf}", file=sys.stderr)
                    return 2
                written, skipped = _gen_from_vcf(vcf, raw_fh, progress_every)
            elif mode == "bed":
                bed = Path(args.bed)
                if not bed.exists():
                    print(f"ERROR: BED not found: {bed}", file=sys.stderr)
                    return 2
                written, skipped = _gen_bed(bed, None, raw_fh, progress_every)
            else:
                written, skipped = _gen_full(args.src, raw_fh, progress_every)

        print(f"Reduced {written} raw rows; sorting + indexing...", file=sys.stderr)
        _sort_file(raw, srt)
        rows = _collapse(srt, fin)
        _index(fin, out)

        meta = {"mode": table_mode, "source": SOURCE, "built_from": built_from, "rows": rows}
        meta_path = out.with_name(out.name + ".meta")
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    finally:
        for t in (raw, srt, fin):
            try:
                t.unlink()
            except OSError:
                pass

    size = out.stat().st_size if out.exists() else 0
    print("-" * 72, file=sys.stderr)
    print(f"Done ({mode}, mode={table_mode}).", file=sys.stderr)
    print(f"  rows written : {rows}", file=sys.stderr)
    if skipped:
        print(f"  skipped      : {skipped} (gnomAD lookup failed after retries)",
              file=sys.stderr)
    print(f"  table        : {out}  ({_human(size)})", file=sys.stderr)
    print(f"  index        : {out}.tbi", file=sys.stderr)
    print(f"  sidecar      : {out}.meta", file=sys.stderr)
    if rows == 0:
        _warn("no rows were written — the table is empty (check network / inputs).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
