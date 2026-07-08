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


def _minimal_alt(ref: str, alt: str) -> str:
    """VEP's minimal ALT: trim the longest shared leading base(s); '-' if emptied."""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    return alt[i:] or "-"


def _allele_match(vep_allele: str, alt: str, ref: str = "") -> bool:
    if vep_allele == alt:
        return True
    if ref:
        return vep_allele == _minimal_alt(ref, alt)
    # No ref available: accept the single-base trim or a fully-trimmed insertion.
    return vep_allele == alt[1:] or vep_allele == "-"


def parse_snpeff(ann: str, alt: str, ref: str = "") -> Optional[dict]:
    first = None
    for entry in ann.split(","):
        f = entry.split("|")
        if len(f) < len(_SNPEFF):
            continue
        if first is None:
            first = f
        if _allele_match(f[0], alt, ref):
            first = f
            break  # first allele-matching entry is most severe for that allele
    if first is None:
        return None
    return {"gene": first[3] or None, "consequence": _first_term(first[1]),
            "hgvs_c": first[9] or None, "hgvs_p": first[10] or None}


def parse_vep(csq: str, alt: str, field_names: list[str], ref: str = "") -> Optional[dict]:
    idx = {name.lower(): i for i, name in enumerate(field_names)}

    def get(f: list[str], name: str) -> Optional[str]:
        i = idx.get(name)
        return f[i] if (i is not None and i < len(f) and f[i]) else None

    def row(f: list[str]) -> dict:
        return {"gene": get(f, "symbol") or get(f, "gene"),
                "consequence": _first_term(get(f, "consequence") or ""),
                "hgvs_c": get(f, "hgvsc"), "hgvs_p": get(f, "hgvsp")}

    entries = [e.split("|") for e in csq.split(",")]
    a = idx.get("allele")

    def matches(f: list[str]) -> bool:
        return a is not None and a < len(f) and _allele_match(f[a], alt, ref)

    # Restrict to rows for THIS ALT first (multiallelic CSQ carries every allele);
    # only if none match (single-allele / no Allele field) consider all rows.
    candidates = [f for f in entries if matches(f)] or entries
    for f in candidates:                      # PICK wins
        if get(f, "pick") == "1":
            return row(f)
    for f in candidates:                      # then CANONICAL
        if get(f, "canonical") == "YES":
            return row(f)
    return row(candidates[0]) if candidates else None


def extract(info: dict[str, str], alt: str,
            csq_format: Optional[list[str]] = None, ref: str = "") -> Optional[dict]:
    """Best consequence/HGVS from ANN, then CSQ, then plain keys."""
    if info.get("ANN"):
        r = parse_snpeff(info["ANN"], alt, ref)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    csq = info.get("CSQ")
    # Only treat CSQ as VEP if it's the structured, pipe-delimited form.
    if csq and "|" in csq and csq_format:
        r = parse_vep(csq, alt, csq_format, ref)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    # Plain keys (synthetic sample / simple pipelines).
    simple = {"gene": info.get("GENE"),
              "consequence": info.get("CSQ") if (info.get("CSQ") and "|" not in info["CSQ"]) else None,
              "hgvs_c": info.get("HGVSC"), "hgvs_p": info.get("HGVSP")}
    return simple if any(simple.values()) else None
