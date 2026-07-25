#!/usr/bin/env python3
"""Build the ClinVar residue index that drives PS1 / PM5.

PS1 (same amino-acid change as a *different* established pathogenic variant) and PM5
(novel missense at a residue where a *different* missense is pathogenic) both need a
lookup from ``(gene, protein_position)`` to the pathogenic amino-acid changes seen
there. ClinVar's per-variant VCF carries no protein change, so this reads the tab-
delimited ``variant_summary.txt.gz`` (whose ``Name`` column is
``NM_…(GENE):c.… (p.Xxx###Yyy)``), keeps GRCh38 pathogenic / likely-pathogenic
missense SNVs with a criteria-based (>=1-star) review, and writes one compact row per
distinct ``(gene, aa_pos, alt_aa)``:

    # gene  aa_pos  ref_aa  alt_aa  stars  genomic_key(chr-pos-ref-alt)  accession

Real, public data (NCBI ClinVar, public domain). The full index is machine-specific
and rebuilt weekly, so it is git-ignored like the other ClinVar store artifacts; a
small committed frozen slice (scripts/freeze_clinvar_residue.py) keeps the offline
demo + tests exercising PS1/PM5 with no download.

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/fetch_clinvar_residue.py
"""
from __future__ import annotations

import gzip
import re
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
OUT = REPO / "data" / "clinvar" / "clinvar_residue.tsv.gz"

# 20 standard amino acids (3-letter). A missense is a substitution between two of these;
# Ter (stop), synonymous (=) and indel/dup/fs 'p.' forms are excluded.
_AA3 = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}
_P_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


def _stars(review: str) -> int:
    r = (review or "").lower().replace("_", " ").strip()
    if "practice guideline" in r:
        return 4
    if "reviewed by expert panel" in r:
        return 3
    if "multiple submitters" in r and "no conflict" in r:
        return 2
    if r.startswith("criteria provided"):
        return 1
    return 0


def _is_plp(sig: str) -> bool:
    s = (sig or "").lower()
    # Pathogenic / Likely pathogenic / Pathogenic-Likely pathogenic; never conflicting/benign.
    return "pathogenic" in s and "conflicting" not in s and "benign" not in s


def main() -> int:
    tmp = REPO / "scratch" / "variant_summary.txt.gz"
    if not tmp.exists():
        tmp.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {URL} ...")
        urllib.request.urlretrieve(URL, tmp)

    # best row per (gene, aa_pos, alt_aa): keep the highest star rating.
    best: dict[tuple, tuple] = {}
    n_rows = n_missense = 0
    with gzip.open(tmp, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        gi, ni, si, ri = idx["GeneSymbol"], idx["Name"], idx["ClinicalSignificance"], idx["ReviewStatus"]
        ai, ci = idx["Assembly"], idx["Chromosome"]
        pvi, refi, alti = idx["PositionVCF"], idx["ReferenceAlleleVCF"], idx["AlternateAlleleVCF"]
        vidi = idx["VariationID"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= vidi or f[ai] != "GRCh38":
                continue
            n_rows += 1
            if not _is_plp(f[si]):
                continue
            stars = _stars(f[ri])
            if stars < 1:
                continue
            m = _P_RE.search(f[ni])
            if not m:
                continue
            ref_aa, pos, alt_aa = m.group(1), m.group(2), m.group(3)
            if ref_aa not in _AA3 or alt_aa not in _AA3 or ref_aa == alt_aa:
                continue
            gene = f[gi]
            if not gene or gene == "-":
                continue
            n_missense += 1
            key = f"{f[ci]}-{f[pvi]}-{f[refi].upper()}-{f[alti].upper()}"
            k = (gene, pos, alt_aa)
            prev = best.get(k)
            if prev is None or stars > prev[0]:
                best[k] = (stars, ref_aa, key, f[vidi])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(OUT, "wt") as w:
        w.write("# ClinVar residue index for PS1/PM5 (GRCh38 P/LP missense, >=1 star).\n")
        w.write("# Columns: gene\taa_pos\tref_aa\talt_aa\tstars\tgenomic_key\taccession\n")
        for (gene, pos, alt_aa), (stars, ref_aa, key, vid) in sorted(best.items()):
            w.write(f"{gene}\t{pos}\t{ref_aa}\t{alt_aa}\t{stars}\t{key}\tVCV{int(vid):09d}\n")
    print(f"scanned {n_rows:,} GRCh38 rows, {n_missense:,} P/LP missense >=1★ "
          f"-> {len(best):,} distinct (gene,residue,alt) -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
