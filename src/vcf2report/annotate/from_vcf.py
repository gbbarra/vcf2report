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


def _pick(x: Optional[str], idx: int) -> Optional[str]:
    """The idx-th comma element of an INFO value; clamps for scalars.

    gnomAD AF/AC/nhomalt and ABraOM AF are Number=A (one value per ALT). After the
    multiallelic split, each Variant carries its ALT index so we read THIS allele's
    value, not allele #1's. A scalar (single value) clamps to itself for any index.
    """
    if x is None:
        return None
    parts = str(x).split(",")
    return parts[idx] if idx < len(parts) else parts[-1]


def _num(x: Optional[str], idx: int = 0):
    p = _pick(x, idx)
    if p is None:
        return None
    try:
        return float(p)
    except ValueError:
        return None


def _int(x: Optional[str], idx: int = 0):
    v = _num(x, idx)
    return int(v) if v is not None else None


def extract(variant: Variant) -> dict:
    """Return the annotation fields found in INFO (only present keys)."""
    info = variant.info or {}
    A = config.INFO_ALIASES
    i = variant.alt_index
    out: dict = {}

    gaf = _first(info, A["gnomad_af"])
    if gaf is not None:
        out["gnomad_af"] = _num(gaf, i)
        out["gnomad_ac"] = _int(_first(info, A["gnomad_ac"]), i)
        out["gnomad_an"] = _int(_first(info, A["gnomad_an"]), i)
        out["gnomad_hom"] = _int(_first(info, A["gnomad_hom"]), i)

    abaf = _first(info, A["abraom_af"])
    if abaf is not None:
        out["abraom_af"] = _num(abaf, i)

    # ClinVar CLNSIG/CLNREVSTAT carry no commas within a value (they use / | _),
    # so per-allele indexing is safe for an un-split multiallelic record.
    sig = _pick(_first(info, A["clinvar_sig"]), i)
    if sig is not None:
        out["clinvar_significance"] = str(sig).replace("_", " ")
        rev = _pick(_first(info, A["clinvar_review"]), i)
        out["clinvar_review_status"] = str(rev).replace("_", " ") if rev else None
        cond = _first(info, A["clinvar_disease"])
        out["clinvar_condition"] = str(cond).replace("_", " ") if cond else None
        out["clinvar_accession"] = _pick(_first(info, A["clinvar_accession"]), i)

    rv = _first(info, A["revel"])
    if rv is not None:
        out["revel"] = _num(rv, i)
    cd = _first(info, A["cadd"])
    if cd is not None:
        out["cadd"] = _num(cd, i)
    return out
