#!/usr/bin/env python3
"""Build the real gnomAD gene-constraint table used for PVS1's LoF-intolerance gate.

Downloads gnomAD's published per-gene LoF constraint metrics (v2.1.1, the standard
pLI / LOEUF table) from the public GCS bucket and exports the two columns the
engine needs — ``gene<TAB>pLI<TAB>LOEUF`` — replacing the tiny hand-curated subset
previously bundled. Real, public data; ships gzipped.

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/fetch_constraint.py
"""
from __future__ import annotations

import gzip
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
URL = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/"
       "constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
OUT = REPO / "data" / "constraint" / "gene_constraint.tsv.gz"


def main() -> int:
    tmp = REPO / "scratch" / "constraint.txt.bgz"
    if not tmp.exists():
        tmp.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL} ...")
        urllib.request.urlretrieve(URL, tmp)

    def _f(v):
        try:
            return f"{float(v):.4g}"
        except (TypeError, ValueError):
            return ""

    OUT.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(tmp, "rt") as fh, gzip.open(OUT, "wt") as w:
        hdr = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(hdr)}
        gi, pi, li = idx["gene"], idx["pLI"], idx["oe_lof_upper"]
        w.write("# gnomAD v2.1.1 LoF constraint (by gene). "
                "Columns: gene\tpLI\tLOEUF(oe_lof_upper)\n")
        seen = set()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= li:
                continue
            gene = p[gi]
            if not gene or gene in seen:
                continue      # keep the first (canonical) row per gene symbol
            seen.add(gene)
            w.write(f"{gene}\t{_f(p[pi])}\t{_f(p[li])}\n")
            n += 1
    print(f"Wrote {n} genes to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
