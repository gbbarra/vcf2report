"""Per-variant quality control.

Applies the documented thresholds (DP, GQ, allele balance, FILTER) and returns
the variants that pass plus the reasons any were dropped, so the report can
show an honest QC funnel.
"""
from __future__ import annotations

from ..config import QC_AB_MAX, QC_AB_MIN, QC_MIN_DP, QC_MIN_GQ
from ..models import Variant


def passes_qc(v: Variant) -> tuple[bool, str]:
    # Carrier gate: the proband must actually carry this ALT allele. zygosity is
    # None for a hom-ref (0/0), a no-call (./.), or a variant present only as a
    # *different* ALT. Non-carriers must never reach candidates/report — not even
    # via the ClinVar P/LP rarity/impact bypass in the filter step.
    if v.zygosity is None:
        return False, "non-carrier (hom-ref / no-call / other allele)"
    if v.filter_status and v.filter_status not in ("PASS", ".", ""):
        return False, f"FILTER={v.filter_status}"
    if v.depth is not None and v.depth < QC_MIN_DP:
        return False, f"DP={v.depth}<{QC_MIN_DP}"
    if v.gq is not None and v.gq < QC_MIN_GQ:
        return False, f"GQ={v.gq}<{QC_MIN_GQ}"
    # The allele-balance window tests a HET expectation (~50%). A half-call ("./1") is
    # recorded as het conservatively — the second allele was never called — so its balance
    # has no expected value to test against, and applying the het window would drop the very
    # carrier the half-call handling exists to keep.
    if v.zygosity == "het" and not v.partial_call and v.allele_balance is not None:
        if not (QC_AB_MIN <= v.allele_balance <= QC_AB_MAX):
            return False, f"AB={v.allele_balance} outside [{QC_AB_MIN},{QC_AB_MAX}]"
    return True, "PASS"


def is_metric_drop(reason: str) -> bool:
    """True when a variant was dropped for a reason worth re-checking against ClinVar.

    The proband carries the allele and the call exists — it was removed on quality, either
    a borderline metric (DP/GQ/AB) or the caller's own FILTER. Neither is reconsidered for
    the candidate list; the point is that such a drop must not be able to delete a *known*
    pathogenic variant from the report in silence (see pipeline._qc_rescue).

    A non-carrier genotype is different in kind: there is no allele to report, so it is
    excluded here.
    """
    return reason.startswith(("DP=", "GQ=", "AB=", "FILTER="))


def apply_qc(variants: list[Variant]) -> tuple[list[Variant], list[tuple[Variant, str]]]:
    kept: list[Variant] = []
    dropped: list[tuple[Variant, str]] = []
    for v in variants:
        ok, reason = passes_qc(v)
        if ok:
            kept.append(v)
        else:
            dropped.append((v, reason))
    return kept, dropped
