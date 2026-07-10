"""Estimate sequencing quality from the parsed VCF.

A VCF holds data only at variant sites, so this reports depth/quality **at those
sites** — a proxy for how well the sample sequenced — plus assay-level signals
(Ti/Tv, het:hom, variant count). It deliberately does NOT claim genome-wide breadth
of coverage: a variants-only VCF cannot provide it (a gVCF or BAM would).
"""
from __future__ import annotations

import statistics
from typing import Iterable, Optional

from ..config import QC_MIN_GQ
from ..models import SeqQuality, Variant

# Transitions: purine<->purine (A<->G) and pyrimidine<->pyrimidine (C<->T).
_TRANSITIONS = ({"A", "G"}, {"C", "T"})


def _pct(hits: int, n: int) -> Optional[float]:
    return round(100 * hits / n, 1) if n else None


def _assay_guess(n: int) -> str:
    """Rough assay class from the variant count (single sample, GRCh38)."""
    if n >= 500_000:
        return "whole-genome-scale"
    if n >= 8_000:
        return "exome / large-panel-scale"
    if n >= 200:
        return "targeted-panel-scale"
    return "small / demo VCF"


def _notes(q: SeqQuality) -> list[str]:
    notes = [
        "Depth is measured only at called variant sites — a proxy for sequencing "
        "quality, not genome-wide breadth of coverage (a variants-only VCF cannot "
        "give breadth; a gVCF or BAM would).",
    ]
    if q.titv is not None:
        notes.append(
            f"Ti/Tv = {q.titv} (expected ~3.0 for exome, ~2.0-2.1 for whole genome; "
            "a much lower value suggests false-positive calls)."
        )
    if q.dp_median is not None and q.dp_median < 20:
        notes.append(
            f"Median depth at variant sites is {q.dp_median}x — below a typical 20-30x "
            "clinical target; interpret low-depth calls with caution."
        )
    return notes


def estimate(variants: Iterable[Variant]) -> SeqQuality:
    """Summarise sequencing quality from all called variants."""
    variants = list(variants)
    q = SeqQuality(n_variants=len(variants), assay_guess=_assay_guess(len(variants)))

    dps = [v.depth for v in variants if v.depth is not None]
    if dps:
        q.n_with_dp = len(dps)
        q.dp_mean = round(statistics.fmean(dps), 1)
        q.dp_median = round(statistics.median(dps), 1)
        q.dp_pct_ge10 = _pct(sum(1 for d in dps if d >= 10), len(dps))
        q.dp_pct_ge20 = _pct(sum(1 for d in dps if d >= 20), len(dps))

    gqs = [v.gq for v in variants if v.gq is not None]
    if gqs:
        q.n_with_gq = len(gqs)
        q.gq_median = round(statistics.median(gqs), 1)
        q.gq_pct_ge20 = _pct(sum(1 for g in gqs if g >= QC_MIN_GQ), len(gqs))

    ti = tv = 0
    for v in variants:
        if len(v.ref) != 1 or len(v.alt) != 1:
            continue
        a, b = v.ref.upper(), v.alt.upper()
        if a not in "ACGT" or b not in "ACGT" or a == b:
            continue
        if {a, b} in _TRANSITIONS:
            ti += 1
        else:
            tv += 1
    q.n_snv = ti + tv
    if tv:
        q.titv = round(ti / tv, 2)

    q.n_het = sum(1 for v in variants if v.zygosity == "het")
    q.n_hom = sum(1 for v in variants if v.zygosity == "hom")
    if q.n_hom:
        q.het_hom_ratio = round(q.n_het / q.n_hom, 2)

    q.notes = _notes(q)
    return q
