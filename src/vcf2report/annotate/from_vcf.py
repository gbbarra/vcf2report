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


def _pick(
    x: Optional[str], idx: int, n_alts: int = 1, per_allele: bool = False
) -> Optional[str]:
    """The idx-th comma element of a multi-value INFO field.

    ``per_allele=True`` marks a VCF ``Number=A`` field — gnomAD AF/AC/AN/nhomalt and
    local cohort AF, one value per ALT. For those, a LONE value at a multiallelic site is not a
    site-wide constant: it means the annotator resolved only one allele, and applying it to
    the others hands ALT #2 ALT #1's frequency (a rare allele inheriting a common one's AF
    can flip BA1/BS1). Return None so the caller falls back to the local snapshot.

    ``per_allele=False`` is for genuinely per-SITE values (gnomAD v4.1's
    ``fafmax_faf95_max``), where a single value legitimately applies to every allele.

    Either way, an array that is >1 long but does NOT cover ``idx`` is a mismatch: None.
    """
    if x is None:
        return None
    parts = str(x).split(",")
    if idx < len(parts):
        return parts[idx]
    if len(parts) != 1:
        return None
    return None if (per_allele and n_alts > 1) else parts[0]


def _num(x: Optional[str], idx: int = 0, n_alts: int = 1, per_allele: bool = False):
    p = _pick(x, idx, n_alts, per_allele)
    if p is None or p == ".":
        return None
    try:
        return float(p)
    except ValueError:
        return None


def _int(x: Optional[str], idx: int = 0, n_alts: int = 1, per_allele: bool = False):
    v = _num(x, idx, n_alts, per_allele)
    return int(v) if v is not None else None


def _multi_num(x: Optional[str], idx: int = 0, n_alts: int = 1):
    """Max numeric from a multi-value predictor field.

    A comma is ambiguous in these fields: dbNSFP writes one value per TRANSCRIPT, while a
    VCF ``Number=A`` field writes one per ALT ALLELE. The only available disambiguator is
    the count — when the comma arity equals the record's ALT count at a multiallelic site,
    the commas are per-allele and only THIS allele's element may be read. Taking the max
    across them instead handed the benign allele of a multiallelic site the pathogenic
    allele's score, inverting the in-silico evidence (PP3 firing on a BP4 allele).
    ``;`` and ``&`` are always per-transcript. ``.`` means missing.
    """
    if x is None:
        return None
    import re

    s = str(x)
    parts = s.split(",")
    if n_alts > 1 and len(parts) == n_alts:
        s = parts[idx] if idx < len(parts) else ""
    vals = []
    for tok in re.split(r"[;,&]", s):
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

    n = variant.n_alts or 1
    gaf = _first(info, A["gnomad_af"])
    if gaf is not None:
        out["gnomad_af"] = _num(gaf, i, n, per_allele=True)
        out["gnomad_ac"] = _int(_first(info, A["gnomad_ac"]), i, n, per_allele=True)
        out["gnomad_an"] = _int(_first(info, A["gnomad_an"]), i, n, per_allele=True)
        out["gnomad_hom"] = _int(_first(info, A["gnomad_hom"]), i, n, per_allele=True)
        # fafmax_faf95_max is a single per-SITE value in gnomAD v4.1 — the maximum filtering
        # AF across the site's alleles. `_benign_af` returns faf95 outright when it is
        # present, ahead of the allele's own AF, so handing a site maximum to a rare ALT
        # called it Benign on another allele's frequency: an allele at AF 1e-5 met BA1 from
        # a site faf95 of 0.29, with a trail that reported PM2 ("absent from population
        # databases") and "0.2900 exceeds 0.05" in the same breath.
        #
        # A filtering AF is derived from one allele's AC/AN, so a site maximum cannot be
        # attributed to a specific allele. Use it only where the site has ONE ALT (where the
        # two are the same thing) or where the array is genuinely per-allele; otherwise leave
        # it unset and let BA1/BS1 fall back to this allele's own AF.
        faf = _first(info, A["gnomad_faf95"])
        if faf is not None:
            per_allele_faf = len(str(faf).split(",")) == n
            if n == 1 or per_allele_faf:
                out["gnomad_faf95"] = _num(faf, i, n, per_allele=per_allele_faf)

    abaf = _first(info, A["local_cohort_af"])
    if abaf is not None:
        out["local_cohort_af"] = _num(abaf, i, n, per_allele=True)

    # ClinVar CLNSIG/CLNREVSTAT contain literal commas (e.g. "Pathogenic,_low_
    # penetrance"), so they cannot be comma-indexed by allele the way the numeric
    # Number=A fields are.
    #
    # At a MULTIALLELIC site that leaves no way to tell which ALT the assertion describes,
    # and applying it to all of them is not a harmless approximation: a 21%-frequency
    # synonymous allele inherited a 2-star Pathogenic assertion from the SNP beside it,
    # bypassed the rarity/impact gate through filter.py's ClinVar P/LP rescue, and was
    # presented to the reader in the report's do-not-dismiss net. The benign direction is
    # symmetric — a novel frameshift inheriting BP6 from a common neighbour.
    #
    # So the site-level value is used only when the record has ONE ALT. Otherwise this
    # returns nothing and annotate falls through to the allele-keyed ClinVar client, which
    # answers for the exact variant or honestly reports no record. Per-allele disambiguation
    # from CLNALLELEID would let INFO be used here too, and remains out of scope.
    sig = _first(info, A["clinvar_sig"]) if n == 1 else None
    if sig is not None:
        out["clinvar_significance"] = str(sig).replace("_", " ")
        rev = _first(info, A["clinvar_review"])
        out["clinvar_review_status"] = str(rev).replace("_", " ") if rev else None
        cond = _first(info, A["clinvar_disease"])
        out["clinvar_condition"] = str(cond).replace("_", " ") if cond else None
        out["clinvar_accession"] = _first(info, A["clinvar_accession"])

    # REVEL/CADD/AlphaMissense are per-transcript multi-values (max = most damaging) UNLESS
    # the comma arity matches the ALT count, in which case they are per-allele — see _multi_num.
    rv = _multi_num(_first(info, A["revel"]), i, n)
    if rv is not None:
        out["revel"] = rv
    cd = _multi_num(_first(info, A["cadd"]), i, n)
    if cd is not None:
        out["cadd"] = cd
    am = _multi_num(_first(info, A["am_pathogenicity"]), i, n)
    if am is not None:
        out["am_pathogenicity"] = am
        amc = (
            _pick(_first(info, A["am_class"]), i)
            if n > 1
            else _first(info, A["am_class"])
        )
        if amc is not None:
            out["am_class"] = str(amc)
    return out
