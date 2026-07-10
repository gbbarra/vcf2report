"""Estimate sequencing quality from the parsed VCF.

A VCF holds data only at variant sites, so this reports depth/quality **at those
sites** — a proxy for how well the sample sequenced — plus assay-level signals
(Ti/Tv, het:hom, variant count). It deliberately does NOT claim genome-wide breadth
of coverage: a variants-only VCF cannot provide it (a gVCF or BAM would).
"""
from __future__ import annotations

import statistics
from typing import Iterable, Optional

from ..config import QC_AB_MAX, QC_AB_MIN, QC_MIN_GQ
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
    if q.n_variants and q.pct_novel is None:
        notes.append(
            "VCF is not dbSNP-annotated (few/no rsIDs in the ID column) — novelty "
            "rate not computed (annotate with dbSNP IDs to enable it)."
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

    # Indel:SNV ratio (indel = a length change; SNV counted above).
    q.n_indel = sum(1 for v in variants if len(v.ref) != len(v.alt))
    if q.n_snv:
        q.indel_snv_ratio = round(q.n_indel / q.n_snv, 3)

    # % multiallelic sites — records are split per ALT, so regroup by site.
    site_alts: dict = {}
    for v in variants:
        k = (v.chrom, v.pos)
        site_alts[k] = max(site_alts.get(k, 1), v.n_alts)
    q.n_sites = len(site_alts)
    q.n_multiallelic_sites = sum(1 for a in site_alts.values() if a > 1)
    q.pct_multiallelic = _pct(q.n_multiallelic_sites, q.n_sites)

    # Novelty vs dbSNP — only meaningful when the VCF is actually dbSNP-annotated
    # (most known variants carry an rsID). A handful of stray rsIDs is not enough.
    q.n_with_rsid = sum(1 for v in variants if (v.variant_id or "").startswith("rs"))
    if variants and q.n_with_rsid >= 0.2 * len(variants):
        q.pct_novel = _pct(len(variants) - q.n_with_rsid, len(variants))

    # Het allele balance in [QC_AB_MIN, QC_AB_MAX] — a contamination / miscall signal.
    het_ab = [v.allele_balance for v in variants
              if v.zygosity == "het" and v.allele_balance is not None]
    if het_ab:
        q.n_het_ab = len(het_ab)
        q.pct_het_ab_balanced = _pct(
            sum(1 for ab in het_ab if QC_AB_MIN <= ab <= QC_AB_MAX), len(het_ab))

    # Fraction of records with FILTER = PASS (or unfiltered).
    if variants:
        q.pct_pass = _pct(
            sum(1 for v in variants if (v.filter_status or "") in ("PASS", ".", "")),
            len(variants))

    q.notes = _notes(q)
    return q
