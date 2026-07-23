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

from .. import config
from ..annotate import extra
from ..annotate.dosage import haploinsufficient
from ..annotate.inheritance import lof_is_disease_mechanism
from ..config import AF_BA1
from ..models import Annotation, CriterionResult, Variant

# Tunable in-silico / frequency cutoffs (documented in the report methods).
REVEL_PATHOGENIC = 0.70
REVEL_BENIGN = 0.15
CADD_PATHOGENIC = 20.0
CADD_BENIGN = 10.0
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
def _is_last_exon(exon: Optional[str]) -> bool:
    """True if the annotation's exon rank 'N/M' places the variant in the last exon
    (N == M). Single-exon transcripts (1/1) also escape NMD, so they count."""
    if not exon or "/" not in exon:
        return False
    n, _, m = exon.partition("/")
    try:
        return int(m) > 0 and int(n) == int(m)
    except ValueError:
        return False


def _pvs1_strength(v: Variant) -> str:
    """Deterministic PVS1 strength via the ClinGen SVI decision tree (Abou Tayoun
    et al. 2018), to the depth a single annotated VCF supports.

    The tree engages ONLY when the annotation carries an exon rank (VEP EXON /
    SnpEff rank); without it PVS1 stays Very Strong. This single gate keeps the
    refinement purely additive on tree-annotated VCFs, so un-annotated inputs — the
    synthetic demos and the frozen panel — are provably unchanged. Real annotators
    always emit an exon rank alongside these consequences (the start codon and any
    premature stop sit in a coding exon), so a genuine clinical variant is graded;
    only degenerate consequence-only records fall back to Very Strong.

      * ``start_lost`` -> Moderate (initiation-codon variant; no downstream evidence)
      * nonsense/frameshift in the LAST exon (NMD-escaping) -> Strong
      * everything else qualifying -> Very Strong (the default)

    The penultimate-exon last-50-nt NMD rule and the >10%-of-protein refinement need
    CDS coordinates we don't carry, so they are intentionally not attempted."""
    if not v.exon:
        return "very_strong"
    cons = v.consequence or ""
    if cons == "start_lost":
        return "moderate"
    if cons in ("stop_gained", "frameshift_variant") and _is_last_exon(v.exon):
        return "strong"
    return "very_strong"


@criterion("PVS1")
def pvs1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Null variant in a gene where LoF is a known disease mechanism"
    # ClinGen SVI's first question is whether LoF is a known disease MECHANISM for the gene.
    # Population constraint (pLI/LOEUF) only measures selection against HETEROZYGOUS LoF, so
    # it is structurally blind to recessive disease: the carrier is healthy, no selection
    # acts, the gene scores "tolerant" — and gating on it rejects the very variants whose
    # mechanism IS loss of function. So an established autosomal-recessive phenotype counts
    # as mechanism evidence alongside constraint. Dominant genes still need the constraint
    # evidence: a dominant phenotype can be gain-of-function or dominant-negative, where a
    # null is NOT the mechanism. Both routes are proxies for ClinGen gene-disease/dosage
    # curation and the report names which one fired.
    ar_mechanism = lof_is_disease_mechanism(v.gene)
    clingen_hi = haploinsufficient(v.gene)
    mechanism = bool(a.gene_lof_intolerant) or clingen_hi or ar_mechanism
    met = bool(v.is_lof and mechanism)
    strength = _pvs1_strength(v) if met else None
    # ClinGen HI=3 is a curated gene-disease-mechanism statement, so it is named first when present;
    # constraint and the recessive-phenotype route are the proxies that fill the gaps it leaves.
    basis = ("gene is ClinGen Haploinsufficiency=3 (curated: LoF causes disease)" if clingen_hi
             else "gene is LoF-intolerant (population constraint)" if a.gene_lof_intolerant
             else "gene has an established autosomal-recessive phenotype (HPO)" if ar_mechanism
             else None)
    cites = []
    if clingen_hi:
        cites.append("ClinGen Dosage Sensitivity (HI=3, local)")
    if a.source.get("gene_lof_intolerant"):
        cites.append(a.source["gene_lof_intolerant"])
    if ar_mechanism and not a.gene_lof_intolerant and not clingen_hi:
        cites.append("HPO gene-to-phenotype inheritance (local)")
    if met and strength != "very_strong":
        reason = (
            f"{v.consequence} is loss-of-function and LoF is a disease mechanism in {v.gene} "
            f"({basis}), but the ClinGen SVI PVS1 tree downgrades it to {strength} "
            f"({'initiation-codon variant' if v.consequence == 'start_lost' else f'NMD-escaping (last exon {v.exon})'})"
        )
    elif met:
        reason = f"{v.consequence} is loss-of-function and LoF is a disease mechanism in {v.gene} ({basis})"
    elif v.is_lof:
        reason = (f"{v.consequence} is a null variant, but LoF is not an established disease "
                  f"mechanism in {v.gene}: not LoF-intolerant by constraint and no known "
                  f"autosomal-recessive phenotype")
    else:
        reason = f"{v.consequence or 'variant'} is not a qualifying null variant"
    return CriterionResult(
        "PVS1", name, "very_strong", applies=True, met=met,
        applied_strength=strength,
        evidence={"consequence": v.consequence, "gene_lof_intolerant": a.gene_lof_intolerant,
                  "lof_mechanism_basis": basis, "gene_moi": config.gene_inheritance(v.gene),
                  "exon": v.exon, "pvs1_strength": strength},
        citation=cites, reasoning=reason,
    )


@criterion("PS1")
def ps1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Same amino-acid change as a DIFFERENT established pathogenic variant"
    # PS1 requires a residue-level cross-match to a *distinct* pathogenic variant
    # (a different nucleotide change giving the same amino-acid change). That needs
    # a residue-indexed ClinVar lookup we don't have, so it is surfaced for expert/
    # model adjudication. The variant's OWN ClinVar assertion is captured by PP5,
    # not PS1 (using it here would double-count reputable-source evidence at Strong).
    return CriterionResult(
        "PS1", name, "strong", applies=True, met=False, adjudicated_by="model",
        confidence="moderate",
        evidence={"hgvs_p": v.hgvs_p},
        reasoning="Requires a residue-level cross-match to a distinct ClinVar "
                  "pathogenic variant — model adjudication (own ClinVar record is PP5)",
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
    ceiling, moi = config.pm2_af_ceiling(v.gene)
    moi_note = f"{v.gene} is {moi}" if moi else "inheritance unknown → strict default"
    rare_global = gaf <= ceiling
    rare_local = baf <= ceiling
    met = rare_global and rare_local
    # Don't quote an ABraOM AF we never checked: None means "not in the local table",
    # not a verified 0.0 (it still doesn't block PM2 — baf defaults to 0.0 above).
    abraom_txt = f"ABraOM AF={a.abraom_af:.6f}" if a.abraom_af is not None else "ABraOM not checked"
    reason = (
        f"gnomAD popmax AF={gaf:.6f}, {abraom_txt} — gnomAD at/under {ceiling:g} ({moi_note})"
        if met else
        f"present above the {ceiling:g} PM2 ceiling ({moi_note}): "
        f"gnomAD AF={gaf:.6f}, {abraom_txt}"
    )
    return CriterionResult(
        "PM2", name, "moderate", applies=True, met=met,
        applied_strength=config.pm2_strength() if met else None,
        evidence={"gnomad_af": a.gnomad_af, "abraom_af": a.abraom_af,
                  "ceiling": ceiling, "moi": moi, "strength_model": config.acmg_model()},
        citation=cites, reasoning=reason,
    )


@criterion("PM3")
def pm3(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PM3", "Detected in trans with a pathogenic variant (recessive)", "moderate",
               "Requires phasing / a second variant — not determinable from this VCF alone")


@criterion("PM4")
def pm4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Protein length change (in-frame indel / stop-loss) in a non-repeat region"
    # Recognise any in-frame indel term (VEP inframe_insertion/deletion, SnpEff
    # disruptive_/conservative_inframe_*, or a generic inframe_indel) + stop-loss.
    c = (v.consequence or "").lower()
    met = "inframe" in c or c == "stop_lost"
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
    name = "Missense in a gene with a low rate of benign missense and where missense is a mechanism"
    # Engine-decided from gnomAD gene-level missense constraint: a gene significantly
    # depleted of missense variation (mis_z >= 3.09, the ClinGen SVI cut) has a low
    # benign-missense rate, so a missense variant there earns PP2 (Supporting). Only
    # missense variants qualify; the metric is absent for many genes -> not met.
    is_missense = (v.consequence or "") == "missense_variant"
    constrained = bool(a.gene_missense_constrained)
    met = is_missense and constrained
    cites = [a.source["gene_constraint"]] if a.source.get("gene_constraint") else []
    mz = a.gene_mis_z
    if met:
        reason = (f"missense in {v.gene}, a missense-constrained gene "
                  f"(gnomAD mis_z={mz:.2f} ≥ {extra.MIS_Z_CONSTRAINED})")
    elif not is_missense:
        reason = f"{v.consequence or 'variant'} is not a missense variant"
    elif mz is not None:
        reason = (f"{v.gene} missense z-score {mz:.2f} is below the "
                  f"{extra.MIS_Z_CONSTRAINED} missense-constraint threshold")
    else:
        reason = f"no gnomAD missense-constraint metric for {v.gene or 'variant'}"
    return CriterionResult(
        "PP2", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        adjudicated_by="engine", confidence="high",
        evidence={"consequence": v.consequence, "gene": v.gene, "mis_z": mz,
                  "mis_z_cutoff": extra.MIS_Z_CONSTRAINED,
                  "missense_constrained": a.gene_missense_constrained},
        citation=cites, reasoning=reason,
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
    # AlphaMissense, calibrated to a variable ACMG strength (ClinGen 2024), takes
    # precedence when available: a strong pathogenic prediction can reach PP3_Strong,
    # which (with PM2) is enough to lift a rare missense from VUS to Likely Pathogenic.
    if a.am_pathogenicity is not None:
        strength = config.am_pp3_strength(a.am_pathogenicity)
        met = strength is not None
        return CriterionResult(
            "PP3", name, "supporting", applies=True, met=met,
            applied_strength=strength if met else None,
            evidence={"am_pathogenicity": a.am_pathogenicity, "am_class": a.am_class,
                      "predictor": "AlphaMissense (ClinGen-calibrated)"},
            citation=[c for c in [a.source.get("alphamissense")] if c],
            confidence="moderate",
            reasoning=(f"AlphaMissense={a.am_pathogenicity:.3f} -> PP3_{strength}"
                       if met else
                       f"AlphaMissense={a.am_pathogenicity:.3f} below the PP3 pathogenic threshold"),
        )
    # Fallback: REVEL/CADD at Supporting strength.
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


@criterion("PP5")
def pp5(v: Variant, a: Annotation) -> CriterionResult:
    name = "Reputable source (ClinVar) classifies the variant as pathogenic"
    sig = (a.clinvar_significance or "").lower()
    is_plp = sig.startswith("pathogenic") or sig.startswith("likely pathogenic")
    # ClinVar review status comes with spaces (E-utilities) or underscores (VCF);
    # normalize, then require a criteria-based (>=1 star) assertion. Use startswith
    # so the 0-star "no assertion criteria provided" (which CONTAINS the substring
    # "criteria provided") is correctly excluded.
    review = (a.clinvar_review_status or "").lower().replace("_", " ").strip()
    reviewed = (review.startswith("criteria provided")
                or "reviewed by expert" in review
                or "practice guideline" in review)
    met = bool(is_plp and reviewed)
    # PP5 was deprecated by the ClinGen SVI; retained here as a transparent, gated
    # SUPPORTING line so ClinVar contributes without over-weighting (vs the old PS1).
    return CriterionResult(
        "PP5", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"clinvar": a.clinvar_significance, "review_status": a.clinvar_review_status},
        citation=[a.clinvar_accession] if (met and a.clinvar_accession) else [],
        confidence="moderate",
        reasoning=(f"ClinVar {a.clinvar_significance} ({a.clinvar_review_status})"
                   if met else "no reviewed ClinVar pathogenic assertion (or 0-star)"),
    )


# ===========================================================================
# Benign criteria
# ===========================================================================
def _benign_af(a: Annotation) -> tuple[float, str]:
    """Allele frequency to test for benign criteria (BA1/BS1), and its basis.

    Prefers gnomAD's filtering AF (faf95, grpmax) — the 95%-CI-bounded value
    ClinGen/Whiffin recommend for frequency-based benign criteria, robust to a
    single small subpopulation inflating a raw popmax. Falls back to the popmax
    AF (max of gnomAD grpmax and ABraOM) when faf95 is not available.
    """
    if a.gnomad_faf95 is not None:
        return a.gnomad_faf95, "gnomAD filtering AF (faf95, grpmax)"
    # None (not 0.0) when NOTHING was looked up -> BA1/BS1 report 'unavailable' rather
    # than a fabricated 0.0 that reads as a checked value (matches PM2's honesty).
    vals = [x for x in (a.gnomad_af, a.abraom_af) if x is not None]
    if not vals:
        return None, "no gnomAD/ABraOM frequency available"
    return max(vals), "gnomAD/ABraOM popmax AF (no faf95 available)"


@criterion("BA1")
def ba1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Allele frequency > 5% in a population database (stand-alone benign)"
    af, basis = _benign_af(a)
    met = af is not None and af >= AF_BA1
    reasoning = (f"{basis} — cannot assess" if af is None
                 else f"{basis} = {af:.4f} exceeds {AF_BA1:g}" if met
                 else f"{basis} = {af:.4f} below {AF_BA1:g}")
    return CriterionResult(
        "BA1", name, "stand_alone", applies=True, met=met,
        applied_strength="stand_alone" if met else None,
        evidence={"af": af, "cutoff": AF_BA1, "basis": basis},
        citation=[c for c in (a.source.get("gnomad"), a.source.get("abraom")) if c],
        reasoning=reasoning,
    )


@criterion("BS1")
def bs1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Allele frequency greater than expected for the disorder"
    af, basis = _benign_af(a)
    cutoff, moi = config.bs1_af_cutoff(v.gene)
    moi_note = f"{v.gene} is {moi}" if moi else "inheritance unknown → default cutoff"
    # BS1 is the "too common for the disorder" band below BA1's stand-alone 5%.
    met = af is not None and cutoff <= af < AF_BA1
    reasoning = (f"{basis} — cannot assess ({moi_note})" if af is None
                 else f"{basis} = {af:.4f} ≥ {cutoff:g} ({moi_note}), below BA1's {AF_BA1:g}" if met
                 else f"{basis} = {af:.4f} under the {cutoff:g} BS1 cutoff ({moi_note})")
    return CriterionResult(
        "BS1", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        evidence={"af": af, "cutoff": cutoff, "moi": moi, "basis": basis},
        citation=[c for c in (a.source.get("gnomad"), a.source.get("abraom")) if c],
        reasoning=reasoning,
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


@criterion("BP1")
def bp1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Missense variant in a gene where primarily truncating variants cause disease"
    # Engine proxy for BP1 (no curated gene list): a missense variant in a gene that is
    # LoF-intolerant (truncating variants are the disease mechanism) yet shows NO missense
    # depletion (gnomAD oe_mis_upper >= 1.0 — missense is tolerated). Both conditions
    # together mean missense is unlikely to be pathogenic in that gene. Deliberately
    # conservative and mutually exclusive with PP2 (a gene cannot be both missense-
    # constrained and missense-tolerant), so no variant earns PP2 and BP1 at once.
    is_missense = (v.consequence or "") == "missense_variant"
    lof_mech = bool(a.gene_lof_intolerant)
    tolerant = bool(a.gene_missense_tolerant) and not a.gene_missense_constrained
    met = is_missense and lof_mech and tolerant
    cites = [a.source["gene_constraint"]] if a.source.get("gene_constraint") else []
    if met:
        reason = (f"missense in {v.gene}: LoF-intolerant gene (truncating is the mechanism) "
                  f"that tolerates missense (gnomAD oe_mis upper ≥ {extra.OE_MIS_TOLERANT})")
    elif not is_missense:
        reason = f"{v.consequence or 'variant'} is not a missense variant"
    elif not lof_mech:
        reason = f"{v.gene or 'gene'} is not LoF-intolerant — truncating-only mechanism not established"
    else:
        reason = f"{v.gene or 'gene'} does not tolerate missense (not depleted-free) — BP1 not supported"
    return CriterionResult(
        "BP1", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        adjudicated_by="engine", confidence="moderate",
        evidence={"consequence": v.consequence, "gene": v.gene,
                  "gene_lof_intolerant": a.gene_lof_intolerant,
                  "oe_mis_upper": a.gene_oe_mis_upper,
                  "oe_mis_cutoff": extra.OE_MIS_TOLERANT,
                  "missense_tolerant": a.gene_missense_tolerant,
                  "missense_constrained": a.gene_missense_constrained},
        citation=cites, reasoning=reason,
    )


@criterion("BP4")
def bp4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Multiple in-silico lines of evidence suggest no impact"
    # AlphaMissense takes precedence when available. Richards Table 5 has no benign
    # "moderate" bucket, so AlphaMissense benign evidence is capped at Supporting.
    if a.am_pathogenicity is not None:
        strength = config.am_bp4_strength(a.am_pathogenicity)
        met = strength is not None
        return CriterionResult(
            "BP4", name, "supporting", applies=True, met=met,
            applied_strength=strength if met else None,
            evidence={"am_pathogenicity": a.am_pathogenicity, "am_class": a.am_class,
                      "predictor": "AlphaMissense (ClinGen-calibrated)"},
            citation=[c for c in [a.source.get("alphamissense")] if c],
            confidence="moderate",
            reasoning=(f"AlphaMissense={a.am_pathogenicity:.3f} -> BP4_{strength}"
                       if met else
                       f"AlphaMissense={a.am_pathogenicity:.3f} above the BP4 benign threshold"),
        )
    # Fallback: REVEL/CADD at Supporting strength.
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


@criterion("BP6")
def bp6(v: Variant, a: Annotation) -> CriterionResult:
    name = "Reputable source (ClinVar) classifies the variant as benign"
    sig = (a.clinvar_significance or "").lower()
    # Mirror of PP5: a criteria-based ClinVar Benign / Likely benign assertion. "conflicting"
    # is deliberately excluded (it starts with neither token) — a conflicted record is not a
    # reputable benign call. The variant's OWN assertion, not a residue cross-match.
    is_blp = sig.startswith("benign") or sig.startswith("likely benign")
    # ClinVar review status arrives space- (E-utilities) or underscore-delimited (VCF); normalize,
    # then require a criteria-based (>=1 star) assertion — same gate as PP5, so 0-star
    # "no assertion criteria provided" is excluded despite containing "criteria provided".
    review = (a.clinvar_review_status or "").lower().replace("_", " ").strip()
    reviewed = (review.startswith("criteria provided")
                or "reviewed by expert" in review
                or "practice guideline" in review)
    met = bool(is_blp and reviewed)
    # BP6 was deprecated by the ClinGen SVI (like PP5); retained here as a transparent, gated
    # SUPPORTING benign line so a reviewed ClinVar benign assertion contributes symmetrically to
    # PP5 without over-weighting. PP5 and BP6 are mutually exclusive (a record is P or B, not both),
    # so no double-count; if benign evidence conflicts with pathogenic criteria the combiner reports VUS.
    return CriterionResult(
        "BP6", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        evidence={"clinvar": a.clinvar_significance, "review_status": a.clinvar_review_status},
        citation=[a.clinvar_accession] if (met and a.clinvar_accession) else [],
        confidence="moderate",
        reasoning=(f"ClinVar {a.clinvar_significance} ({a.clinvar_review_status})"
                   if met else "no reviewed ClinVar benign assertion (or 0-star)"),
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
