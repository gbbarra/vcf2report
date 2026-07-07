"""Variant tiering: turn thousands of calls into a ranked candidate shortlist.

Each funnel step is recorded so the shortlist is explainable in the report.
Runs on already-annotated variants (needs population AF + ClinVar + HPO match).
"""
from __future__ import annotations

from dataclasses import dataclass

from ..config import AF_RECESSIVE_MAX
from ..models import Annotation, Variant

# Consequences we keep as clinically relevant (coding / canonical splice).
IMPACTFUL = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
    "missense_variant", "inframe_insertion", "inframe_deletion",
}


@dataclass
class FilterFunnel:
    total: int = 0
    after_rarity: int = 0
    after_impact: int = 0
    candidates: int = 0
    notes: list[str] = None  # type: ignore
    # Variants dropped as common in ABraOM despite being rare/absent in gnomAD —
    # i.e. spurious candidates a gnomAD-only pipeline would have kept. This is the
    # concrete, per-run evidence for the Brazilian-frequency differentiator.
    abraom_filtered: list[str] = None  # type: ignore

    def __post_init__(self):
        if self.notes is None:
            self.notes = []
        if self.abraom_filtered is None:
            self.abraom_filtered = []


def _is_clinvar_plp(a: Annotation) -> bool:
    sig = (a.clinvar_significance or "").lower()
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


def _is_rare(a: Annotation, max_af: float) -> bool:
    af = max(a.gnomad_af or 0.0, a.abraom_af or 0.0)
    return af <= max_af


def filter_variants(
    annotated: list[tuple[Variant, Annotation]],
    max_af: float = AF_RECESSIVE_MAX,
) -> tuple[list[tuple[Variant, Annotation]], FilterFunnel]:
    """Return (ranked candidates, funnel). ClinVar P/LP bypass rarity/impact."""
    funnel = FilterFunnel(total=len(annotated))

    # Step 1 — rarity (ClinVar P/LP always retained regardless of AF).
    rare = [(v, a) for v, a in annotated if _is_rare(a, max_af) or _is_clinvar_plp(a)]
    funnel.after_rarity = len(rare)

    # Record variants dropped by ABraOM that a gnomAD-only filter would keep:
    # rare/absent in gnomAD but common (> cutoff) in the Brazilian cohort.
    for v, a in annotated:
        gnomad_rare = (a.gnomad_af or 0.0) <= max_af
        abraom_common = (a.abraom_af or 0.0) > max_af
        if gnomad_rare and abraom_common and not _is_clinvar_plp(a):
            funnel.abraom_filtered.append(
                f"{v.gene or v.key} ({v.hgvs_p or v.key}): gnomAD AF={a.gnomad_af or 0.0:.6f} "
                f"but ABraOM AF={a.abraom_af:.4f} — common in Brazilians, dropped"
            )

    # Step 2 — impact (coding/splice; ClinVar P/LP retained regardless).
    impactful = [
        (v, a) for v, a in rare
        if (v.consequence in IMPACTFUL) or _is_clinvar_plp(a)
    ]
    funnel.after_impact = len(impactful)

    # Step 3 — phenotype ranking: on-phenotype variants and ClinVar P/LP float up.
    def rank_key(pair: tuple[Variant, Annotation]) -> tuple:
        v, a = pair
        return (
            0 if _is_clinvar_plp(a) else 1,          # ClinVar P/LP first
            -(a.hpo_match_score or 0.0),             # higher phenotype match first
            max(a.gnomad_af or 0.0, a.abraom_af or 0.0),  # rarer first
        )

    candidates = sorted(impactful, key=rank_key)
    funnel.candidates = len(candidates)
    funnel.notes.append(
        f"{funnel.total} variants -> {funnel.after_rarity} rare -> "
        f"{funnel.after_impact} coding/splice -> {funnel.candidates} candidates"
    )
    return candidates, funnel
