"""Parse consequence/HGVS from real annotator output in a VCF's INFO.

Supports SnpEff ``ANN`` and Ensembl VEP ``CSQ`` (whose subfield order is read
from the VCF header). Falls back to plain ``GENE``/``CSQ``/``HGVSC``/``HGVSP``
keys (as the bundled synthetic sample uses). Returns a dict with gene,
consequence, hgvs_c, hgvs_p — or None if nothing usable is found.
"""
from __future__ import annotations

import re
from typing import Optional

# SnpEff ANN subfield order (fixed by the SnpEff spec).
_SNPEFF = ["allele", "annotation", "impact", "gene", "gene_id", "feature_type",
           "feature_id", "biotype", "rank", "hgvs_c", "hgvs_p"]


def parse_csq_format(header_lines: list[str]) -> Optional[list[str]]:
    """Extract the VEP CSQ subfield names from the ##INFO CSQ header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ") and "Format:" in line:
            fmt = line.split("Format:")[1].strip().rstrip('">').strip()
            return [f.strip() for f in fmt.split("|")]
    return None


def _first_term(consequence: str) -> Optional[str]:
    # VEP/SnpEff join multiple consequences with "&"; the first is most severe.
    return (consequence.split("&")[0].strip() or None) if consequence else None


def _allele_match(vep_allele: str, alt: str) -> bool:
    # SNVs match directly; VEP trims the shared leading base for indels.
    return vep_allele == alt or vep_allele == alt[1:] or vep_allele == "-"


def parse_snpeff(ann: str, alt: str) -> Optional[dict]:
    first = None
    for entry in ann.split(","):
        f = entry.split("|")
        if len(f) < len(_SNPEFF):
            continue
        if first is None:
            first = f
        if _allele_match(f[0], alt):
            first = f
            break  # first allele-matching entry is most severe for that allele
    if first is None:
        return None
    return {"gene": first[3] or None, "consequence": _first_term(first[1]),
            "hgvs_c": first[9] or None, "hgvs_p": first[10] or None}


def parse_vep(csq: str, alt: str, field_names: list[str]) -> Optional[dict]:
    idx = {name.lower(): i for i, name in enumerate(field_names)}

    def get(f: list[str], name: str) -> Optional[str]:
        i = idx.get(name)
        return f[i] if (i is not None and i < len(f) and f[i]) else None

    def row(f: list[str]) -> dict:
        return {"gene": get(f, "symbol") or get(f, "gene"),
                "consequence": _first_term(get(f, "consequence") or ""),
                "hgvs_c": get(f, "hgvsc"), "hgvs_p": get(f, "hgvsp")}

    entries = [e.split("|") for e in csq.split(",")]
    canonical = first = allele_hit = None
    for f in entries:
        if first is None:
            first = f
        if get(f, "pick") == "1":
            return row(f)
        if canonical is None and get(f, "canonical") == "YES":
            canonical = f
        a = idx.get("allele")
        if allele_hit is None and a is not None and a < len(f) and _allele_match(f[a], alt):
            allele_hit = f
    chosen = canonical or allele_hit or first
    return row(chosen) if chosen else None


def extract(info: dict[str, str], alt: str,
            csq_format: Optional[list[str]] = None) -> Optional[dict]:
    """Best consequence/HGVS from ANN, then CSQ, then plain keys."""
    if info.get("ANN"):
        r = parse_snpeff(info["ANN"], alt)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    csq = info.get("CSQ")
    # Only treat CSQ as VEP if it's the structured, pipe-delimited form.
    if csq and "|" in csq and csq_format:
        r = parse_vep(csq, alt, csq_format)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    # Plain keys (synthetic sample / simple pipelines).
    simple = {"gene": info.get("GENE"),
              "consequence": info.get("CSQ") if (info.get("CSQ") and "|" not in info["CSQ"]) else None,
              "hgvs_c": info.get("HGVSC"), "hgvs_p": info.get("HGVSP")}
    return simple if any(simple.values()) else None
