#!/usr/bin/env python3
"""Slice the full local ClinVar residue index down to a small committed frozen slice.

``scripts/fetch_clinvar_residue.py`` builds the genome-wide index (git-ignored, ~MBs).
This keeps only the rows for a curated set of genes — the worked-example and e2e demo
genes by default — and writes ``data/clinvar/clinvar_residue_frozen.tsv.gz``, which IS
committed so PS1/PM5 fire out-of-the-box on the shipped examples with no download.

    python scripts/freeze_clinvar_residue.py                 # default demo gene set
    python scripts/freeze_clinvar_residue.py --genes BRCA1 TP53 SCN1A
    python scripts/freeze_clinvar_residue.py --in data/clinvar/clinvar_residue.tsv.gz
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_IN = REPO / "data" / "clinvar" / "clinvar_residue.tsv.gz"
OUT = REPO / "data" / "clinvar" / "clinvar_residue_frozen.tsv.gz"

# Genes exercised by the committed examples (data/example) and the e2e regression suite.
DEMO_GENES = [
    "SCN1A",
    "RB1",
    "KCNQ2",
    "APC",
    "SCN2A",
    "STK11",
    "STXBP1",
    "WT1",
    "SLC2A1",
    "FBN1",
    "NIPBL",
    "PIGA",
    "BBS2",
    "TGFBR1",
    "SPINT2",
    "RBSN",
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Freeze a small ClinVar residue slice for the demo."
    )
    ap.add_argument(
        "--in", dest="src", default=str(DEFAULT_IN), help="full residue index (.tsv.gz)"
    )
    ap.add_argument("--out", default=str(OUT), help="frozen slice output (.tsv.gz)")
    ap.add_argument(
        "--genes", nargs="+", default=DEMO_GENES, help="gene symbols to keep"
    )
    ap.add_argument(
        "--all",
        action="store_true",
        help="keep EVERY gene (ship the residue index genome-wide, ~1 MB) instead of "
        "the demo subset — PS1/PM5 then fire for all genes out-of-the-box",
    )
    args = ap.parse_args(argv)

    src = Path(args.src)
    if not src.exists():
        raise SystemExit(f"missing {src} — run scripts/fetch_clinvar_residue.py first")
    keep = None if args.all else set(args.genes)

    kept = genes_seen = 0
    seen = set()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    label = "genome-wide" if args.all else "FROZEN demo slice"
    with gzip.open(src, "rt") as fh, gzip.open(args.out, "wt") as w:
        w.write(f"# ClinVar residue index (PS1/PM5) — {label}.\n")
        w.write(
            "# Columns: gene\taa_pos\tref_aa\talt_aa\tstars\tgenomic_key\taccession\n"
        )
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            gene = line.split("\t", 1)[0]
            if keep is None or gene in keep:
                w.write(line)
                kept += 1
                seen.add(gene)
    print(f"wrote {kept:,} rows for {len(seen)} genes -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
