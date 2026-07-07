"""ACMG/AMP criterion evaluators (Richards et al., Genet Med 2015).

Each function evaluates ONE criterion and returns a :class:`CriterionResult`
carrying the concrete evidence values it used, the source citation, and a
one-line reason. This is what makes the classification auditable: no criterion
flips on without a trail you can read.

Two classes of criteria:

* **engine** — deterministic lookups/thresholds (population AF, LoF mechanics,
  in-silico cutoffs). Reproducible and unit-tested.
* **model** — genuine clinical judgment (hotspot, functional-domain, functional
  studies). We surface the evidence and leave the on/off decision tagged
  ``adjudicated_by="model"`` so Claude adjudicates it transparently rather than
  the engine pretending it's a fact.

Criteria requiring data we don't have from a single proband VCF (trio,
segregation, phasing) are returned with ``applies=False`` and an explicit
reason, so the output is honest instead of silently incomplete.
"""
from __future__ import annotations

from typing import Callable

from ..config import AF_BA1, AF_DOMINANT_MAX
from ..models import Annotation, CriterionResult, Variant

# Tunable in-silico / frequency cutoffs (documented in the report methods).
REVEL_PATHOGENIC = 0.70
REVEL_BENIGN = 0.15
CADD_PATHOGENIC = 20.0
CADD_BENIGN = 10.0
PM2_RARE_AF = 1e-4          # "absent or ultra-rare" ceiling
BS1_COMMON_AF = 0.01        # more common than expected for a rare disorder
BS2_HOM_MIN = 2             # healthy homozygotes incompatible with severe disease
HPO_PP4_MIN = 0.60          # phenotype-match score to support PP4

CriterionFn = Callable[[Variant, Annotation], CriterionResult]

_REGISTRY: dict[str, CriterionFn] = {}


def criterion(code: str) -> Callable[[CriterionFn], CriterionFn]:
    def deco(fn: CriterionFn) -> CriterionFn:
        _REGISTRY[code] = fn
        return fn
    return deco


def all_criteria() -> dict[str, CriterionFn]:
    return dict(_REGISTRY)


def _na(code: str, name: str, strength: str, reason: str) -> CriterionResult:
    """A criterion that does not apply to a single-proband VCF."""
    return CriterionResult(
        code=code, name=name, default_strength=strength, applies=False,
        met=False, reasoning=reason, confidence="high", adjudicated_by="engine",
    )


# ===========================================================================
# Pathogenic criteria
# ===========================================================================
@criterion("PVS1")
def pvs1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Null variant in a gene where LoF is a known disease mechanism"
    met = bool(v.is_lof and a.gene_lof_intolerant)
    cites = []
    if a.source.get("gene_lof_intolerant"):
        cites.append(a.source["gene_lof_intolerant"])
    reason = (
        f"{v.consequence} is loss-of-function and {v.gene} is LoF-intolerant"
        if met else
        f"{v.consequence or 'variant'} is not a qualifying null variant in a LoF-intolerant gene"
    )
    return CriterionResult(
        "PVS1", name, "very_strong", applies=True, met=met,
        applied_strength="very_strong" if met else None,
        evidence={"consequence": v.consequence, "gene_lof_intolerant": a.gene_lof_intolerant},
        citation=cites, reasoning=reason,
    )


@criterion("PS1")
def ps1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Same amino-acid change as an established pathogenic variant"
    # Deterministic proxy: ClinVar reports this exact protein change as Pathogenic.
    met = bool(v.hgvs_p and (a.clinvar_significance or "").lower().startswith("pathogenic"))
    cites = [a.clinvar_accession] if (met and a.clinvar_accession) else []
    return CriterionResult(
        "PS1", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        evidence={"hgvs_p": v.hgvs_p, "clinvar": a.clinvar_significance},
        citation=[c for c in cites if c],
        reasoning=(
            f"ClinVar records {v.hgvs_p} as {a.clinvar_significance}"
            if met else "no established pathogenic record for this amino-acid change"
        ),
    )


@criterion("PS2")
def ps2(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PS2", "De novo (confirmed) in a patient", "strong",
               "Requires parental (trio) data — not available from a single proband VCF")


@criterion("PS3")
def ps3(v: Variant, a: Annotation) -> CriterionResult:
    return CriterionResult(
        "PS3", "Well-established functional studies show a damaging effect",
        "strong", applies=True, met=False, adjudicated_by="model",
        confidence="low",
        reasoning="Requires literature review of functional assays — left for expert/model adjudication",
    )


@criterion("PS4")
def ps4(v: Variant, a: Annotation) -> CriterionResult:
    return CriterionResult(
        "PS4", "Prevalence in affected significantly increased vs controls",
        "strong", applies=True, met=False, adjudicated_by="model", confidence="low",
        evidence={"gnomad_af": a.gnomad_af, "abraom_af": a.abraom_af},
        reasoning="Needs case-control data; population absence alone is captured by PM2",
    )


@criterion("PM1")
def pm1(v: Variant, a: Annotation) -> CriterionResult:
    return CriterionResult(
        "PM1", "Located in a mutational hotspot / critical functional domain",
        "moderate", applies=True, met=False, adjudicated_by="model", confidence="moderate",
        evidence={"consequence": v.consequence, "hgvs_p": v.hgvs_p},
        reasoning="Domain/hotspot membership requires curated annotation — model adjudication",
    )


@criterion("PM2")
def pm2(v: Variant, a: Annotation) -> CriterionResult:
    """Absent / ultra-rare in population databases — gnomAD AND ABraOM.

    Checking ABraOM (Brazilian SABE cohort) alongside gnomAD is the key local
    value-add: a variant absent from gnomAD but common in admixed Brazilians
    must NOT earn PM2, which prevents a real class of local misclassifications.
    """
    name = "Absent or ultra-rare in population databases (gnomAD + ABraOM)"
    # gnomAD AF None means 'frequency unavailable' (lookup failed), NOT absence —
    # PM2 must not fire because we cannot assert the variant is rare.
    gnomad_unknown = a.gnomad_af is None
    baf = a.abraom_af if a.abraom_af is not None else 0.0
    cites = [c for c in (a.source.get("gnomad"), a.source.get("abraom")) if c]
    if gnomad_unknown:
        return CriterionResult(
            "PM2", name, "moderate", applies=True, met=False, adjudicated_by="engine",
            confidence="low",
            evidence={"gnomad_af": None, "abraom_af": a.abraom_af},
            citation=cites,
            reasoning="gnomAD frequency unavailable — cannot assert population absence",
        )
    gaf = a.gnomad_af
    rare_global = gaf <= PM2_RARE_AF
    rare_local = baf <= PM2_RARE_AF
    met = rare_global and rare_local
    reason = (
        f"gnomAD popmax AF={gaf:.6f}, ABraOM AF={baf:.6f} — both at/under {PM2_RARE_AF:g}"
        if met else
        f"present in a population DB (gnomAD AF={gaf:.6f}, ABraOM AF={baf:.6f})"
    )
    return CriterionResult(
        "PM2", name, "moderate", applies=True, met=met,
        applied_strength="moderate" if met else None,
        evidence={"gnomad_af": a.gnomad_af, "abraom_af": a.abraom_af},
        citation=cites, reasoning=reason,
    )


@criterion("PM3")
def pm3(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PM3", "Detected in trans with a pathogenic variant (recessive)", "moderate",
               "Requires phasing / a second variant — not determinable from this VCF alone")


@criterion("PM4")
def pm4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Protein length change (in-frame indel / stop-loss) in a non-repeat region"
    met = (v.consequence or "") in {"inframe_insertion", "inframe_deletion", "stop_lost"}
    return CriterionResult(
        "PM4", name, "moderate", applies=True, met=met,
        applied_strength="moderate" if met else None,
        evidence={"consequence": v.consequence},
        reasoning=(f"{v.consequence} alters protein length" if met
                   else "no protein-length-changing consequence"),
    )


@criterion("PM5")
def pm5(v: Variant, a: Annotation) -> CriterionResult:
    return CriterionResult(
        "PM5", "Novel missense at a residue where a different pathogenic missense is known",
        "moderate", applies=True, met=False, adjudicated_by="model", confidence="moderate",
        evidence={"hgvs_p": v.hgvs_p},
        reasoning="Requires residue-level ClinVar cross-check — model adjudication",
    )


@criterion("PM6")
def pm6(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PM6", "Assumed de novo (parentage not confirmed)", "moderate",
               "Requires parental data — not available from a single proband VCF")


@criterion("PP2")
def pp2(v: Variant, a: Annotation) -> CriterionResult:
    return CriterionResult(
        "PP2", "Missense in a gene with a low rate of benign missense and where missense is a mechanism",
        "supporting", applies=True, met=False, adjudicated_by="model", confidence="moderate",
        evidence={"consequence": v.consequence, "gene": v.gene},
        reasoning="Gene-level missense constraint requires curated metric — model adjudication",
    )


def _insilico_direction(a: Annotation) -> Optional[str]:
    """'pathogenic' | 'benign' | 'conflicting' | None from REVEL/CADD.

    PP3 and BP4 are mutually exclusive: if predictors disagree (one deleterious,
    one benign) neither fires, so a variant can never earn both a pathogenic- and
    a benign-supporting line from the same in-silico evidence.
    """
    patho = (a.revel is not None and a.revel >= REVEL_PATHOGENIC) or \
            (a.cadd_phred is not None and a.cadd_phred >= CADD_PATHOGENIC)
    benign = (a.revel is not None and a.revel <= REVEL_BENIGN) or \
             (a.cadd_phred is not None and a.cadd_phred <= CADD_BENIGN)
    if patho and benign:
        return "conflicting"
    if patho:
        return "pathogenic"
    if benign:
        return "benign"
    return None


@criterion("PP3")
def pp3(v: Variant, a: Annotation) -> CriterionResult:
    name = "Multiple in-silico lines of evidence support a deleterious effect"
    direction = _insilico_direction(a)
    met = direction == "pathogenic"
    return CriterionResult(
        "PP3", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"revel": a.revel, "cadd_phred": a.cadd_phred,
                  "revel_cutoff": REVEL_PATHOGENIC, "cadd_cutoff": CADD_PATHOGENIC},
        citation=[c for c in [a.source.get("insilico")] if c],
        reasoning=(f"REVEL={a.revel}, CADD={a.cadd_phred} above deleterious cutoffs"
                   if met else ("in-silico predictors conflict — neither PP3 nor BP4 applied"
                                if _insilico_direction(a) == "conflicting"
                                else "in-silico predictors below deleterious cutoffs / unavailable")),
    )


@criterion("PP4")
def pp4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Patient phenotype highly specific for the gene (HPO match)"
    score = a.hpo_match_score if a.hpo_match_score is not None else 0.0
    met = score >= HPO_PP4_MIN
    return CriterionResult(
        "PP4", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"hpo_match_score": a.hpo_match_score, "matched_terms": a.hpo_matched_terms,
                  "cutoff": HPO_PP4_MIN},
        citation=[c for c in [a.source.get("hpo")] if c],
        reasoning=(f"phenotype match {score:.2f} (terms: {', '.join(a.hpo_matched_terms) or 'n/a'})"
                   if met else f"phenotype match {score:.2f} below {HPO_PP4_MIN}"),
    )


# ===========================================================================
# Benign criteria
# ===========================================================================
@criterion("BA1")
def ba1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Allele frequency > 5% in a population database (stand-alone benign)"
    af = max(a.gnomad_af or 0.0, a.abraom_af or 0.0)
    met = af >= AF_BA1
    return CriterionResult(
        "BA1", name, "stand_alone", applies=True, met=met,
        applied_strength="stand_alone" if met else None,
        evidence={"max_af": af, "cutoff": AF_BA1},
        citation=[c for c in (a.source.get("gnomad"), a.source.get("abraom")) if c],
        reasoning=(f"max population AF={af:.4f} exceeds {AF_BA1:g}"
                   if met else f"max population AF={af:.4f} below {AF_BA1:g}"),
    )


@criterion("BS1")
def bs1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Allele frequency greater than expected for the disorder"
    af = max(a.gnomad_af or 0.0, a.abraom_af or 0.0)
    met = AF_DOMINANT_MAX < af < AF_BA1 and af >= BS1_COMMON_AF
    return CriterionResult(
        "BS1", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        evidence={"max_af": af, "cutoff": BS1_COMMON_AF},
        reasoning=(f"AF={af:.4f} is common relative to a rare-disorder expectation"
                   if met else f"AF={af:.4f} not in the BS1 window"),
    )


@criterion("BS2")
def bs2(v: Variant, a: Annotation) -> CriterionResult:
    name = "Observed in healthy adult homozygotes (incompatible with severe early-onset disease)"
    homs = a.gnomad_homozygotes or 0
    met = homs >= BS2_HOM_MIN
    return CriterionResult(
        "BS2", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        evidence={"gnomad_homozygotes": homs, "cutoff": BS2_HOM_MIN},
        citation=[c for c in [a.source.get("gnomad")] if c],
        reasoning=(f"{homs} homozygotes in gnomAD" if met
                   else f"{homs} homozygotes (below {BS2_HOM_MIN})"),
    )


@criterion("BP4")
def bp4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Multiple in-silico lines of evidence suggest no impact"
    direction = _insilico_direction(a)
    met = direction == "benign"
    return CriterionResult(
        "BP4", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"revel": a.revel, "cadd_phred": a.cadd_phred,
                  "revel_cutoff": REVEL_BENIGN, "cadd_cutoff": CADD_BENIGN},
        reasoning=(f"REVEL={a.revel}, CADD={a.cadd_phred} below benign cutoffs"
                   if met else ("in-silico predictors conflict — neither PP3 nor BP4 applied"
                                if direction == "conflicting"
                                else "in-silico predictors not benign / unavailable")),
    )


@criterion("BP7")
def bp7(v: Variant, a: Annotation) -> CriterionResult:
    name = "Synonymous variant with no predicted splice impact"
    met = (v.consequence or "") == "synonymous_variant"
    return CriterionResult(
        "BP7", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"consequence": v.consequence},
        reasoning=("synonymous with no splice prediction" if met
                   else "not a synonymous variant"),
    )
