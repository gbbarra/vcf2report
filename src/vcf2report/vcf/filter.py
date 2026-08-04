"""Variant tiering: turn thousands of calls into a ranked candidate shortlist.

Each funnel step is recorded so the shortlist is explainable in the report.
Runs on already-annotated variants (needs population AF + ClinVar + HPO match).
"""

from __future__ import annotations

from dataclasses import dataclass

from ..config import AF_RECESSIVE_MAX
from ..models import Annotation, Variant

# Consequences we keep as clinically relevant (coding / canonical splice). In-frame indels
# are spelled differently across annotators (VEP: inframe_insertion/deletion; SnpEff:
# disruptive_/conservative_inframe_*), so recognise any term containing "inframe" rather than
# an exact allowlist — otherwise a real in-frame indel is silently dropped before classification.
IMPACTFUL = {
    "stop_gained",
    "frameshift_variant",
    "splice_donor_variant",
    "splice_acceptor_variant",
    "start_lost",
    "stop_lost",
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    # Whole-transcript / whole-exon loss. These are the MOST damaging consequences an
    # annotator emits, and omitting them meant annparse translated SnpEff's EXON_DELETED
    # into "transcript_ablation" only for this filter to discard it one step later.
    "transcript_ablation",
    "exon_loss_variant",
    # VEP's catch-all for a protein-length/sequence change it cannot type more precisely.
    "protein_altering_variant",
}


# Splice-adjacent terms deliberately kept OUT of IMPACTFUL. They are NOT the canonical
# +/-1,2 sites (those are splice_donor_variant / splice_acceptor_variant, which ARE
# impactful) — they sit 3-8 bases into the intron or in the last 3 exonic bases. Excluding
# them is a real sensitivity limit, so the funnel counts and reports it rather than leaving
# the reader to assume the shortlist covered everything near a splice junction.
#
# Two facts make the exclusion defensible, both measured on the committed exomes:
#   * a KNOWN pathogenic splice-region variant already reaches the shortlist through the
#     ClinVar P/LP bypass in filter_variants, so nothing with evidence is lost;
#   * admitting them adds ~28% more candidates, every one of which lands at VUS — no
#     splice predictor (SpliceAI/MaxEntScan) is wired in, so the engine has no evidence
#     with which to raise them.
# Wire a splice predictor and this trade-off should be revisited.
NEAR_SPLICE = {
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}


def is_inframe_indel(consequence) -> bool:
    """True for any protein-length-changing in-frame indel term (VEP, SnpEff, or generic)."""
    return "inframe" in (consequence or "").lower()


def is_impactful(consequence) -> bool:
    return consequence in IMPACTFUL or is_inframe_indel(consequence)


@dataclass
class FilterFunnel:
    total: int = 0
    after_rarity: int = 0
    after_impact: int = 0
    candidates: int = 0
    notes: list[str] = None  # type: ignore
    # Variants dropped as common in local cohort despite being rare/absent in gnomAD —
    # i.e. spurious candidates a gnomAD-only pipeline would have kept. This is the
    # concrete, per-run evidence that the local cohort earned its place in the funnel.
    local_cohort_filtered: list[str] = None  # type: ignore
    # Rare variants set aside at the impact step because they are splice-ADJACENT rather
    # than canonical splice or coding (see NEAR_SPLICE). Counted so the report states the
    # limit instead of implying the shortlist covered every splice-relevant variant.
    near_splice_excluded: int = 0

    def __post_init__(self):
        if self.notes is None:
            self.notes = []
        if self.local_cohort_filtered is None:
            self.local_cohort_filtered = []


def _is_clinvar_plp(a: Annotation) -> bool:
    # Underscores normalized first: raw ClinVar CLNSIG is underscore-delimited
    # ("Likely_pathogenic"). This test guards the rarity/impact RESCUE, so failing it
    # silently drops a known pathogenic variant from the candidate list entirely.
    sig = (a.clinvar_significance or "").lower().replace("_", " ")
    return sig.startswith("pathogenic") or sig.startswith("likely pathogenic")


def _is_rare(a: Annotation, max_af: float) -> bool:
    af = max(a.gnomad_af or 0.0, a.local_cohort_af or 0.0)
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

    # Record variants dropped by local cohort that a gnomAD-only filter would keep:
    # rare/absent in gnomAD but common (> cutoff) in the operator's cohort.
    for v, a in annotated:
        gnomad_rare = (a.gnomad_af or 0.0) <= max_af
        local_cohort_common = (a.local_cohort_af or 0.0) > max_af
        if gnomad_rare and local_cohort_common and not _is_clinvar_plp(a):
            funnel.local_cohort_filtered.append(
                f"{v.gene or v.key} ({v.hgvs_p or v.key}): gnomAD AF={a.gnomad_af or 0.0:.6f} "
                f"but local cohort AF={a.local_cohort_af:.4f} — common locally, dropped"
            )

    # Step 2 — impact (coding/splice; ClinVar P/LP retained regardless).
    impactful = [
        (v, a) for v, a in rare if is_impactful(v.consequence) or _is_clinvar_plp(a)
    ]
    funnel.after_impact = len(impactful)
    funnel.near_splice_excluded = sum(
        1 for v, a in rare if v.consequence in NEAR_SPLICE and not _is_clinvar_plp(a)
    )

    # Step 3 — phenotype ranking: on-phenotype variants and ClinVar P/LP float up.
    def rank_key(pair: tuple[Variant, Annotation]) -> tuple:
        v, a = pair
        return (
            0 if _is_clinvar_plp(a) else 1,  # ClinVar P/LP first
            -(a.hpo_match_score or 0.0),  # higher phenotype match first
            max(a.gnomad_af or 0.0, a.local_cohort_af or 0.0),  # rarer first
        )

    candidates = sorted(impactful, key=rank_key)
    funnel.candidates = len(candidates)
    funnel.notes.append(
        f"{funnel.total} variants -> {funnel.after_rarity} rare -> "
        f"{funnel.after_impact} coding/splice -> {funnel.candidates} candidates"
    )
    return candidates, funnel
