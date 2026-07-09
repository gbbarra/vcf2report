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
    """The idx-th comma element of a Number=A INFO value.

    gnomAD AF/AC/nhomalt and ABraOM AF are one value per ALT. A genuine scalar
    (single value) is used for any allele; but if the array length is >1 and does
    NOT cover ``idx`` (annotator/ALT-count mismatch), return None rather than
    silently broadcasting allele #1's value onto another allele.
    """
    if x is None:
        return None
    parts = str(x).split(",")
    if idx < len(parts):
        return parts[idx]
    return parts[0] if len(parts) == 1 else None


def _num(x: Optional[str], idx: int = 0):
    p = _pick(x, idx)
    if p is None or p == ".":
        return None
    try:
        return float(p)
    except ValueError:
        return None


def _int(x: Optional[str], idx: int = 0):
    v = _num(x, idx)
    return int(v) if v is not None else None


def _multi_num(x: Optional[str]):
    """Max numeric from a multi-value predictor field (dbNSFP REVEL/CADD are
    per-transcript, separated by ',', ';' or '&'; '.' means missing)."""
    if x is None:
        return None
    import re
    vals = []
    for tok in re.split(r"[;,&]", str(x)):
        tok = tok.strip()
        if tok and tok != ".":
            try:
                vals.append(float(tok))
            except ValueError:
                pass
    return max(vals) if vals else None


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
        # fafmax_faf95_max is a single per-site value in gnomAD v4.1 (not per-allele);
        # _num tolerates both a scalar and an allele-indexed array.
        faf = _first(info, A["gnomad_faf95"])
        if faf is not None:
            out["gnomad_faf95"] = _num(faf, i)

    abaf = _first(info, A["abraom_af"])
    if abaf is not None:
        out["abraom_af"] = _num(abaf, i)

    # ClinVar CLNSIG/CLNREVSTAT contain literal commas (e.g. "Pathogenic,_low_
    # penetrance"), so do NOT comma-index by allele — take the whole value. (Real
    # per-allele ClinVar disambiguation would key on CLNALLELEID, out of scope.)
    sig = _first(info, A["clinvar_sig"])
    if sig is not None:
        out["clinvar_significance"] = str(sig).replace("_", " ")
        rev = _first(info, A["clinvar_review"])
        out["clinvar_review_status"] = str(rev).replace("_", " ") if rev else None
        cond = _first(info, A["clinvar_disease"])
        out["clinvar_condition"] = str(cond).replace("_", " ") if cond else None
        out["clinvar_accession"] = _first(info, A["clinvar_accession"])

    # REVEL/CADD are per-transcript multi-values, not per-allele: aggregate (max).
    rv = _multi_num(_first(info, A["revel"]))
    if rv is not None:
        out["revel"] = rv
    cd = _multi_num(_first(info, A["cadd"]))
    if cd is not None:
        out["cadd"] = _num(cd, i)
    return out
