"""Read annotations already present in a VCF's INFO (SnpEff/VEP + vcfanno).

When a real exome is annotated upstream (the recommended production flow), the
gnomAD AF, ClinVar, and in-silico scores are in INFO. Reading them here means the
whole pipeline runs offline with zero per-variant DB lookups — the fast path for
a real exome. Field names are resolved via ``config.INFO_ALIASES``.
"""
from __future__ import annotations

from typing import Optional

from .. import config
from ..models import Variant


def _first(info: dict, keys: list[str]) -> Optional[str]:
    for k in keys:
        v = info.get(k)
        if v not in (None, ".", ""):
            return v
    return None


def _num(x: Optional[str]):
    if x is None:
        return None
    try:
        return float(str(x).split(",")[0])  # first allele if array-valued
    except ValueError:
        return None


def _int(x: Optional[str]):
    v = _num(x)
    return int(v) if v is not None else None


def extract(variant: Variant) -> dict:
    """Return the annotation fields found in INFO (only present keys)."""
    info = variant.info or {}
    A = config.INFO_ALIASES
    out: dict = {}

    gaf = _first(info, A["gnomad_af"])
    if gaf is not None:
        out["gnomad_af"] = _num(gaf)
        out["gnomad_ac"] = _int(_first(info, A["gnomad_ac"]))
        out["gnomad_an"] = _int(_first(info, A["gnomad_an"]))
        out["gnomad_hom"] = _int(_first(info, A["gnomad_hom"]))

    abaf = _first(info, A["abraom_af"])
    if abaf is not None:
        out["abraom_af"] = _num(abaf)

    sig = _first(info, A["clinvar_sig"])
    if sig is not None:
        out["clinvar_significance"] = str(sig).replace("_", " ")
        rev = _first(info, A["clinvar_review"])
        out["clinvar_review_status"] = str(rev).replace("_", " ") if rev else None
        cond = _first(info, A["clinvar_disease"])
        out["clinvar_condition"] = str(cond).replace("_", " ") if cond else None
        out["clinvar_accession"] = _first(info, A["clinvar_accession"])

    rv = _first(info, A["revel"])
    if rv is not None:
        out["revel"] = _num(rv)
    cd = _first(info, A["cadd"])
    if cd is not None:
        out["cadd"] = _num(cd)
    return out
