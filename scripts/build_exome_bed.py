#!/usr/bin/env python3
"""Build a kit-agnostic hg38 exome BED — MANE (or all protein-coding) exons ±pad — from GENCODE.

Unlike a vendor capture panel (licensed, kit-specific), this is an open, reproducible
target for slicing gnomAD to the exome (scripts/build_gnomad_parquet.py --bed), so
vcf2report's local frequency store covers *any* exome sample, not one lab's panel.
Default = MANE (the NCBI/EMBL-EBI clinical-standard one-transcript-per-gene set, public
domain), read via GENCODE's ``tag=MANE_Select`` / ``tag=MANE_Plus_Clinical`` on exons —
so we get MANE with correct chr-prefixed coords from a single GENCODE download.

    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_exome_bed.py                  # MANE
    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_exome_bed.py --select protein_coding
    python3 scripts/build_exome_bed.py --gff gencode.v46.annotation.gff3.gz        # local

Output: data/gnomad/exome_hg38.bed (chrom, start, end; 0-based half-open; sorted+merged).
"""
from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               "release_46/gencode.v46.basic.annotation.gff3.gz")
_CHROMS = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}


def _open(gff: str):
    op = gzip.open if str(gff).endswith(".gz") else open
    return op(gff, "rt")


def exon_intervals(gff: str, select: str):
    """Yield (chrom, start0, end) for each exon on a standard contig, filtered by
    ``select``: 'mane' (MANE_Select + MANE_Plus_Clinical), 'protein_coding', or 'all'."""
    with _open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in _CHROMS:
                continue
            a = f[8]                                       # tags are a comma-list in tag=...
            if select == "mane" and "MANE_Select" not in a and "MANE_Plus_Clinical" not in a:
                continue
            if select == "protein_coding" and "gene_type=protein_coding" not in a:
                continue
            yield f[0], int(f[3]) - 1, int(f[4])          # GFF3 is 1-based inclusive


def merge(intervals, pad: int):
    """Sort, pad ±pad, and merge overlapping/adjacent intervals per chrom."""
    by_chrom: dict[str, list] = {}
    for chrom, s, e in intervals:
        by_chrom.setdefault(chrom, []).append((max(0, s - pad), e + pad))
    out = []
    for chrom in sorted(by_chrom, key=lambda c: (len(c), c)):
        spans = sorted(by_chrom[chrom])
        cs, ce = spans[0]
        for s, e in spans[1:]:
            if s <= ce:                                    # overlap/adjacent -> extend
                ce = max(ce, e)
            else:
                out.append((chrom, cs, ce)); cs, ce = s, e
        out.append((chrom, cs, ce))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gff", help="local GENCODE GFF3(.gz); else download the default.")
    ap.add_argument("--pad", type=int, default=50, help="bp flanking each exon (default 50).")
    ap.add_argument("--select", choices=["mane", "protein_coding", "all"], default="mane",
                    help="which exons (default: mane — clinical-standard, kit-agnostic).")
    ap.add_argument("--out", default=str(Path(config.GNOMAD_LOCAL_TABIX).parent / "exome_hg38.bed"))
    args = ap.parse_args(argv)

    gff = args.gff
    if not gff:
        if config.offline():
            print("ERROR: needs network to fetch GENCODE (or pass --gff). "
                  "Set VCF2REPORT_ALLOW_NETWORK=1.", file=sys.stderr)
            return 2
        gff = str(Path(args.out).parent / "gencode.basic.gff3.gz")
        Path(gff).parent.mkdir(parents=True, exist_ok=True)
        print(f"Downloading {GENCODE_URL} ...", file=sys.stderr)
        urllib.request.urlretrieve(GENCODE_URL, gff)

    print(f"Extracting + merging {args.select} exons ...", file=sys.stderr)
    merged = merge(exon_intervals(gff, args.select), args.pad)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    span = sum(e - s for _c, s, e in merged)
    with open(out, "w") as fh:
        for chrom, s, e in merged:
            fh.write(f"{chrom}\t{s}\t{e}\n")
    print(f"Wrote {out}: {len(merged):,} regions, {span/1e6:.1f} Mb "
          f"(pad ±{args.pad}, {args.select})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
