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
    if v.zygosity == "het" and v.allele_balance is not None:
        if not (QC_AB_MIN <= v.allele_balance <= QC_AB_MAX):
            return False, f"AB={v.allele_balance} outside [{QC_AB_MIN},{QC_AB_MAX}]"
    return True, "PASS"


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
