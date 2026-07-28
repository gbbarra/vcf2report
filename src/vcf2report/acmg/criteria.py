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

from typing import Callable, Optional

from .. import config
from ..annotate import clinvar_residue as extra_residue
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

# Shared audit-trail wording for the three residue-index criteria (PS1/PM5/PM1). Written once
# because these are user-facing report text that must read identically across the three — and
# because the first embeds a SCRIPT PATH, which would otherwise have to be found in three literals
# if the script is ever renamed.
_RESIDUE_UNAVAILABLE = ("ClinVar residue index unavailable — {code} not assessed "
                        "(build it: scripts/fetch_clinvar_residue.py)")
_OWN_PLP_WITHHELD = ("variant's own ClinVar assertion is pathogenic — captured by PP5; "
                     "{code} withheld to avoid double-counting the same ClinVar evidence")

# The criteria decided from GENE- or RESIDUE-level missense context rather than from the variant's
# own record: gnomAD missense constraint (PP2/BP1) and the ClinVar residue index (PS1/PM5/PM1).
# Defined here, beside the criteria themselves, so adding one is a single edit in a single file —
# the report layer imports this rather than keeping its own copy, which went stale immediately.
MISSENSE_CONTEXT_CODES = ("PS1", "PM1", "PM5", "PP2", "BP1")

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


def _judgment(code: str, name: str, strength: str, reason: str,
              evidence: dict | None = None) -> CriterionResult:
    """A criterion needing genuine clinical/literature judgement.

    The evidence is surfaced and the decision tagged ``adjudicated_by="model"``, defaulting to
    not-met, so a human or model decides it transparently instead of the engine pretending it is a
    fact. One helper so the contract lives in one place rather than being restated per criterion.
    """
    return CriterionResult(
        code=code, name=name, default_strength=strength, applies=True, met=False,
        adjudicated_by="model", confidence="low", evidence=evidence or {}, reasoning=reason,
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
    if a.source.get("gene_constraint"):
        cites.append(a.source["gene_constraint"])
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
        # Only cite when the criterion fired: an unmet PVS1 was still printing "ClinGen Dosage
        # Sensitivity (HI=3)" in the Source column of a row whose reasoning is about consequence
        # type, asserting a curated lookup as the basis of a decision that never used it.
        citation=cites if met else [], reasoning=reason,
    )


def _clinvar_reviewed(a: Annotation) -> bool:
    """Is the ClinVar assertion criteria-based (>=1 star)?

    Review status arrives with spaces (E-utilities) or underscores (VCF); normalize, then use
    ``startswith`` so the 0-star "no assertion criteria provided" — which CONTAINS the substring
    "criteria provided" — is correctly excluded.
    """
    review = (a.clinvar_review_status or "").lower().replace("_", " ").strip()
    return (review.startswith("criteria provided")
            or "reviewed by expert" in review
            or "practice guideline" in review)


def _clinvar_says(a: Annotation, *prefixes: str) -> bool:
    """Does the variant's own ClinVar significance start with one of ``prefixes``?"""
    return (a.clinvar_significance or "").lower().startswith(prefixes)


def _own_clinvar_plp(a: Annotation) -> bool:
    """True when PP5 fires — i.e. the variant's OWN ClinVar record is a REVIEWED P/LP assertion.

    Such a variant is already covered by PP5, so the residue-index criteria PS1/PM5 stand down:
    firing them too would double-count the same ClinVar evidence (the concern behind the SVI's PP5
    deprecation). PS1/PM5 are thus reserved for their real purpose — elevating missense ClinVar has
    NOT directly classified.

    This is deliberately the exact condition ``pp5`` tests, sharing both halves rather than
    restating them: a 0-star "Pathogenic" record does NOT fire PP5, so withholding PS1/PM5 on it
    would strip a legitimate residue cross-match while claiming a PP5 that never fired.
    """
    return _clinvar_says(a, "pathogenic", "likely pathogenic") and _clinvar_reviewed(a)


@criterion("PS1")
def ps1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Same amino-acid change as a DIFFERENT established pathogenic variant"
    # Engine-decided from the ClinVar residue index: a P/LP (>=1-star) missense with the
    # SAME amino-acid change at a DIFFERENT genomic locus. The variant's OWN ClinVar
    # assertion is PP5, not PS1 (the index excludes the query's own key AND a variant already
    # called P/LP is withheld here), so the two never double-count reputable-source evidence.
    m = a.clinvar_ps1
    own_plp = _own_clinvar_plp(a)
    met = m is not None and not own_plp
    if met:
        acc = m.get("accession")
        reason = (f"same amino-acid change ({m.get('ref_aa','')}"
                  f"{'→' + m['alt_aa'] if m.get('alt_aa') else ''}) as an established ClinVar "
                  f"pathogenic variant {acc or ''} ({m.get('stars')}★) at a different locus")
    elif m is not None and own_plp:
        reason = _OWN_PLP_WITHHELD.format(code="PS1")
    elif not a.clinvar_residue_available:
        reason = _RESIDUE_UNAVAILABLE.format(code="PS1")
    elif not v.hgvs_p:
        reason = "no protein change (not a missense) — PS1 not applicable"
    else:
        reason = "no distinct ClinVar pathogenic variant with the same amino-acid change"
    return CriterionResult(
        "PS1", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        adjudicated_by="engine", confidence="high" if a.clinvar_residue_available else "low",
        evidence={"hgvs_p": v.hgvs_p, "ps1_match": m, "own_clinvar_plp": own_plp},
        citation=[m["accession"]] if (met and m.get("accession")) else [],
        reasoning=reason,
    )


@criterion("PS2")
def ps2(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PS2", "De novo (confirmed) in a patient", "strong",
               "Requires parental (trio) data — not available from a single proband VCF")


@criterion("PS3")
def ps3(v: Variant, a: Annotation) -> CriterionResult:
    return _judgment(
        "PS3", "Well-established functional studies show a damaging effect", "strong",
        "Requires literature review of functional assays — left for expert/model adjudication")


@criterion("PS4")
def ps4(v: Variant, a: Annotation) -> CriterionResult:
    return _judgment(
        "PS4", "Prevalence in affected significantly increased vs controls", "strong",
        "Needs case-control data; population absence alone is captured by PM2",
        evidence={"gnomad_af": a.gnomad_af, "abraom_af": a.abraom_af})


def _pm1_signals(v: Variant, a: Annotation) -> dict:
    """Every input to the PM1 decision, plus the decision itself, computed ONCE.

    ``pm1`` needs the individual signals to explain *why* it did or did not fire, and ``pp2`` needs
    only the verdict — so both read this rather than each restating the guards. Keeping the verdict
    next to the signals it is derived from is what stops ``met`` and the reasoning string drifting
    apart when a guard or threshold changes.

    A missense whose NEIGHBOURHOOD is dense with pathogenic missense — and materially denser than
    the gene's own baseline — sits in a mutational hot spot. Withheld when the variant's own residue
    already carries pathogenic evidence (that is PS1/PM5) or when the gene tolerates missense.
    """
    hot = a.clinvar_hotspot or {}
    s = {
        "hot": hot,
        "n_residues": hot.get("n_residues", 0),
        "n_changes": hot.get("n_changes", 0),
        "enrichment": hot.get("enrichment", 0.0),
        "is_missense": (v.consequence or "") == "missense_variant",
        # PM1 stands down only when PS1/PM5 ACTUALLY FIRE, not merely when the index holds a
        # same-residue match. Both are themselves withheld when the variant's own ClinVar record is
        # a reviewed P/LP, so testing raw presence made all three stand down together — and ADDING
        # a pathogenic ClinVar assertion LOWERED the tier (LP -> VUS), which no reading defends.
        # PM1's evidence is the NEIGHBOURING residues (hotspot() excludes the query residue), so it
        # is independent of PP5 and cannot double-count it.
        "same_residue": ((a.clinvar_ps1 is not None or a.clinvar_pm5 is not None)
                         and not _own_clinvar_plp(a)),
        # None (gene absent from the constraint table, or the caller never populated it) is NOT the
        # same as False. `bool(None)` made "unknown" behave exactly like "proven constrained", so
        # PM1's stand-in for ACMG's "without benign variation" silently vanished — the criterion
        # fired with one of its three stated guards inoperative and said nothing about it.
        "tolerant": a.gene_missense_tolerant,
        "tolerance_known": a.gene_missense_tolerant is not None,
    }
    s["dense"] = s["n_residues"] >= extra_residue.HOTSPOT_MIN_RESIDUES
    s["enriched"] = s["enrichment"] >= extra_residue.HOTSPOT_MIN_ENRICHMENT
    # Require the tolerance guard to be ANSWERABLE, not merely non-True: firing PM1 while the
    # constraint metric is unknown means asserting a hotspot with the "without benign variation"
    # half of the criterion never evaluated.
    s["fires"] = bool(s["is_missense"] and s["dense"] and s["enriched"]
                      and not s["same_residue"]
                      and s["tolerance_known"] and not s["tolerant"])
    return s


def _pm1_fires(v: Variant, a: Annotation) -> bool:
    """Whether PM1 applies — the verdict alone, for PP2's stand-down check."""
    return _pm1_signals(v, a)["fires"]


@criterion("PM1")
def pm1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Located in a mutational hotspot / critical functional domain"
    # Engine-decided from the ClinVar residue index: a missense whose NEIGHBOURHOOD (+/- a few
    # residues, the query's own residue excluded) is dense with pathogenic missense is, empirically,
    # in a mutational hot spot. Three guards keep it conservative and non-overlapping:
    #   * PS1/PM5 take precedence — they already carry the SAME-residue evidence, and ACMG warns
    #     against counting PM1 together with PM5. PM1 therefore covers only the residue that is
    #     itself unremarkable but sits inside a pathogenic-dense region.
    #   * "without benign variation" (the second half of the ACMG wording) has no residue-level
    #     benign index here, so it is approximated at gene level: a gene that TOLERATES missense
    #     (BP1's flag) is excluded.
    #   * a high residue count is required (not a single neighbour), so isolated pairs don't fire.
    s = _pm1_signals(v, a)
    hot, n_res, n_chg, enrich = s["hot"], s["n_residues"], s["n_changes"], s["enrichment"]
    is_missense, same_residue_evidence = s["is_missense"], s["same_residue"]
    tolerant, dense, enriched = s["tolerant"], s["dense"], s["enriched"]
    met = s["fires"]
    if met:
        reason = (f"{n_res} distinct pathogenic-missense residues within ±{hot.get('window')} aa "
                  f"({n_chg} pathogenic changes), {enrich}× denser than {v.gene}'s own baseline — "
                  f"mutational hotspot by ClinVar density")
    elif not is_missense:
        reason = f"{v.consequence or 'variant'} is not a missense variant"
    elif not a.clinvar_residue_available:
        reason = _RESIDUE_UNAVAILABLE.format(code="PM1")
    elif same_residue_evidence:
        reason = ("pathogenic missense at this exact residue — carried by PS1/PM5; "
                  "PM1 withheld so the same residue evidence is not counted twice")
    elif tolerant:
        reason = f"{v.gene} tolerates missense (gnomAD obs/exp) — not a constrained hotspot"
    elif not s["tolerance_known"]:
        reason = (f"no gnomAD missense-constraint metric for {v.gene or 'this gene'} — PM1 needs it "
                  f"to stand in for ACMG's 'without benign variation' clause, so it is not assessed")
    elif dense and not enriched:
        reason = (f"{n_res} pathogenic-missense residues nearby, but only {enrich}× {v.gene}'s "
                  f"baseline density (needs {extra_residue.HOTSPOT_MIN_ENRICHMENT}×) — a "
                  f"well-catalogued gene, not a local hotspot")
    else:
        reason = (f"only {n_res} pathogenic-missense residue(s) within ±{hot.get('window', 0)} aa "
                  f"(needs {extra_residue.HOTSPOT_MIN_RESIDUES})")
    return CriterionResult(
        "PM1", name, "moderate", applies=True, met=met,
        applied_strength="moderate" if met else None,
        adjudicated_by="engine", confidence="moderate",
        evidence={"consequence": v.consequence, "hgvs_p": v.hgvs_p,
                  "hotspot_residues": n_res, "hotspot_changes": n_chg,
                  "enrichment": enrich, "gene_baseline": hot.get("gene_baseline"),
                  "window": hot.get("window"), "cutoff": extra_residue.HOTSPOT_MIN_RESIDUES,
                  "enrichment_cutoff": extra_residue.HOTSPOT_MIN_ENRICHMENT},
        citation=["ClinVar residue index (local)"] if met else [],
        reasoning=reason,
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


def _pm5_strength(m: Optional[dict]) -> Optional[str]:
    """Grade PM5 by how well-established the residue is (ClinGen SVI allows PM5 at a variable
    strength rather than a flat Moderate):

    Strength needs BOTH breadth (how many distinct changes at the position) and quality (how well
    reviewed they are) — two single-submitter 1* records at one residue are suggestive, not Strong:

      * **Strong**   — >=2 DISTINCT pathogenic amino-acid changes AND a well-reviewed (>=2*) record:
                       the position, not one substitution, is established as intolerant.
      * **Moderate** — the default: >=2 changes at 1*, or one change with a >=2* record.
      * **Supporting** — a single other pathogenic change resting on one 1* submitter.
    """
    if not m:
        return None
    n_other, stars = (m.get("n_other") or 1), (m.get("stars") or 0)
    if n_other >= 2 and stars >= 2:
        return "strong"
    if n_other >= 2 or stars >= 2:
        return "moderate"
    return "supporting"


@criterion("PM5")
def pm5(v: Variant, a: Annotation) -> CriterionResult:
    name = "Novel missense at a residue where a different pathogenic missense is known"
    # Engine-decided from the ClinVar residue index: a P/LP (>=1-star) missense with a
    # DIFFERENT amino-acid change at the SAME residue, applied only when the query's exact
    # change is not itself established (that case is PS1/PP5). PS1 and PM5 are therefore
    # mutually exclusive — the index sets pm5 to None whenever the same change is known.
    m = a.clinvar_pm5
    own_plp = _own_clinvar_plp(a)
    met = m is not None and a.clinvar_ps1 is None and not own_plp
    strength = _pm5_strength(m) if met else None
    if met:
        acc = m.get("accession")
        n_other = m.get("n_other") or 1
        extra_note = (f"; {n_other} distinct pathogenic changes at this residue → PM5_{strength}"
                      if strength != "moderate" else "")
        reason = (f"a different pathogenic missense at the same residue is established in "
                  f"ClinVar (→{m['alt_aa']}, {acc or ''}, {m.get('stars')}★); this change is novel"
                  f"{extra_note}")
    # own_plp is tested BEFORE the PS1 branch: when the variant's own record is a reviewed P/LP,
    # PS1 is itself withheld, so crediting it would point the reviewer at a criterion that reads
    # "—". Ordering the other way made PM5 say "captured by PS1" while PS1.met was False.
    elif own_plp:
        reason = _OWN_PLP_WITHHELD.format(code="PM5")
    elif a.clinvar_ps1 is not None:
        reason = "same amino-acid change is itself established — captured by PS1, not PM5"
    elif not a.clinvar_residue_available:
        reason = _RESIDUE_UNAVAILABLE.format(code="PM5")
    elif not v.hgvs_p:
        reason = "no protein change (not a missense) — PM5 not applicable"
    else:
        reason = "no other ClinVar pathogenic missense at this residue"
    return CriterionResult(
        "PM5", name, "moderate", applies=True, met=met,
        applied_strength=strength,
        adjudicated_by="engine", confidence="high" if a.clinvar_residue_available else "low",
        evidence={"hgvs_p": v.hgvs_p, "pm5_match": m, "own_clinvar_plp": own_plp,
                  "pm5_strength": strength},
        citation=[m["accession"]] if (met and m.get("accession")) else [],
        reasoning=reason,
    )


@criterion("PM6")
def pm6(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PM6", "Assumed de novo (parentage not confirmed)", "moderate",
               "Requires parental data — not available from a single proband VCF")


@criterion("PP1")
def pp1(v: Variant, a: Annotation) -> CriterionResult:
    return _na("PP1", "Co-segregation with disease in multiple affected family members",
               "supporting",
               "Requires genotypes for affected relatives — not available from a single proband VCF")


@criterion("PP2")
def pp2(v: Variant, a: Annotation) -> CriterionResult:
    name = "Missense in a gene with a low rate of benign missense and where missense is a mechanism"
    # Engine-decided from gnomAD gene-level missense constraint: a gene significantly
    # depleted of missense variation (mis_z >= 3.09, the ClinGen SVI cut) has a low
    # benign-missense rate, so a missense variant there earns PP2 (Supporting). Only
    # missense variants qualify; the metric is absent for many genes -> not met.
    is_missense = (v.consequence or "") == "missense_variant"
    constrained = bool(a.gene_missense_constrained)
    # PP2 and PM1 make the same claim — "missense matters here" — at two granularities (whole gene
    # vs. this region). ClinGen VCEP specifications generally direct that they not both be applied,
    # so the more specific, stronger evidence wins: when PM1 fires, PP2 stands down.
    pm1_fires = _pm1_fires(v, a)
    met = is_missense and constrained and not pm1_fires
    cites = [a.source["gene_constraint"]] if a.source.get("gene_constraint") else []
    mz = a.gene_mis_z
    if met:
        # gene_missense_constrained and gene_mis_z are independent optional fields. The shipped
        # annotator derives the flag FROM the score so they cannot disagree, but a hand-built or
        # deserialised Annotation can set the flag alone — and formatting None with :.2f raised
        # TypeError out of evaluate_criteria, aborting the whole classification rather than
        # degrading one criterion's wording.
        score = f"gnomAD mis_z={mz:.2f} ≥ {extra.MIS_Z_CONSTRAINED}" if mz is not None \
            else "flagged missense-constrained by the annotation (mis_z not supplied)"
        reason = f"missense in {v.gene}, a missense-constrained gene ({score})"
    elif pm1_fires:
        reason = ("regional hotspot evidence applies (PM1, Moderate) — PP2 stands down so the "
                  "same missense-intolerance signal is not counted at two granularities")
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


def _insilico_names(a: Annotation) -> str:
    """Name only the predictors that actually returned a value.

    The reason strings printed "REVEL=None, CADD=32.0" whenever one predictor was absent — and
    REVEL is missense-only, so that is the norm for every LoF/synonymous variant. `None` reads as
    a rendering failure rather than "not computed", and it hid how many lines of evidence the
    criterion really rested on.
    """
    parts = [f"{n}={x}" for n, x in (("REVEL", a.revel), ("CADD", a.cadd_phred)) if x is not None]
    return ", ".join(parts) if parts else "no in-silico predictor available"


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
        reasoning=(f"{_insilico_names(a)} above deleterious cutoffs"
                   if met else ("in-silico predictors conflict — neither PP3 nor BP4 applied"
                                if _insilico_direction(a) == "conflicting"
                                else "in-silico predictors below deleterious cutoffs / unavailable")),
    )


@criterion("PP4")
def pp4(v: Variant, a: Annotation) -> CriterionResult:
    name = "Patient phenotype highly specific for the gene (HPO match)"
    # A score of None means no comparison happened — no HPO terms were supplied, or the gene has no
    # HPO curation. Rendering that as "phenotype match 0.00 below 0.6" reads as a MEASUREMENT that
    # argues against the variant, when nothing was ever compared. `run_pipeline(..., hpo_terms=[])`
    # is a supported call, so a whole report could carry that fabricated line on every candidate.
    score = a.hpo_match_score
    met = score is not None and score >= HPO_PP4_MIN
    return CriterionResult(
        "PP4", name, "supporting", applies=True, met=met,
        applied_strength="supporting" if met else None,
        confidence="high" if score is not None else "low",
        evidence={"hpo_match_score": score, "matched_terms": a.hpo_matched_terms,
                  "cutoff": HPO_PP4_MIN},
        citation=[c for c in [a.source.get("hpo")] if c] if score is not None else [],
        reasoning=("no phenotype comparison — no HPO terms supplied, or the gene has no HPO "
                   "annotation" if score is None
                   else f"phenotype match {score:.2f} (terms: {', '.join(a.hpo_matched_terms) or 'n/a'})"
                   if met else f"phenotype match {score:.2f} below {HPO_PP4_MIN}"),
    )


@criterion("PP5")
def pp5(v: Variant, a: Annotation) -> CriterionResult:
    name = "Reputable source (ClinVar) classifies the variant as pathogenic"
    met = _own_clinvar_plp(a)   # a REVIEWED (>=1-star) P/LP assertion — see the helper
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
    single small subpopulation inflating a raw popmax.

    But it takes the MAXIMUM of that and ABraOM, never faf95 alone. gnomAD has almost no
    admixed-Brazilian representation, and all three gnomAD backends report an absent variant as
    ``faf95 = 0.0`` rather than None — so returning faf95 the moment it exists discarded the
    Brazilian frequency entirely. A variant carried by 20% of Brazilians and absent from gnomAD
    earned neither BA1 nor BS1, while the trail asserted "= 0.0000 below 0.05". Installing the
    local gnomAD store made classification strictly WORSE, which is the opposite of the intent.
    PM2 already reads ABraOM this way; the benign criteria were the inconsistent ones.
    """
    # None (not 0.0) when NOTHING was looked up -> BA1/BS1 report 'unavailable' rather
    # than a fabricated 0.0 that reads as a checked value (matches PM2's honesty).
    faf, braz = a.gnomad_faf95, a.abraom_af
    if faf is not None and braz is not None:
        if braz > faf:
            return braz, "ABraOM (SABE) AF — above the gnomAD filtering AF"
        return faf, "gnomAD filtering AF (faf95, grpmax)"
    if faf is not None:
        return faf, "gnomAD filtering AF (faf95, grpmax)"
    vals = [x for x in (a.gnomad_af, braz) if x is not None]
    if not vals:
        return None, "no gnomAD/ABraOM frequency available"
    return max(vals), "gnomAD/ABraOM popmax AF (no faf95 available)"


@criterion("BA1")
def ba1(v: Variant, a: Annotation) -> CriterionResult:
    name = "Allele frequency > 5% in a population database (stand-alone benign)"
    af, basis = _benign_af(a)
    # Strictly greater: Richards 2015 BA1 is "allele frequency is >5%", and this criterion's own
    # name and met-wording both say "exceeds". At exactly 5.00% `>=` awarded stand-alone Benign —
    # the strongest possible under-call — on a boundary the guideline does not authorise.
    met = af is not None and af > AF_BA1
    reasoning = (f"{basis} — cannot assess" if af is None
                 else f"{basis} = {af:.4f} exceeds {AF_BA1:g}" if met
                 else f"{basis} = {af:.4f} at or below {AF_BA1:g}")
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
    # BS1 is the "too common for the disorder" BAND, bounded above by BA1's stand-alone 5%.
    met = af is not None and cutoff <= af <= AF_BA1
    # Three outcomes, not two: below the band, inside it, or above it (where BA1 takes over). The
    # ladder previously had no arm for "above", so an 8% allele was described as "under the 0.01
    # cutoff" on the line directly beneath BA1 reporting it exceeded 0.05.
    reasoning = (f"{basis} — cannot assess ({moi_note})" if af is None
                 else f"{basis} = {af:.4f} ≥ {cutoff:g} ({moi_note}), at or below BA1's {AF_BA1:g}" if met
                 else f"{basis} = {af:.4f} exceeds BA1's {AF_BA1:g} — superseded by BA1 "
                      f"(stand-alone benign), which subsumes BS1" if af > AF_BA1
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
    # None means the homozygote count was never looked up (no store, no record, or a build
    # mismatch that skipped gnomAD) — NOT an observed zero. `or 0` collapsed the two, so the trail
    # asserted "0 homozygotes in gnomAD" at confidence=high, citing a source that never supplied
    # the field. PM2/BA1/BS1 already distinguish these; BS2 was the outlier.
    homs = a.gnomad_homozygotes
    met = homs is not None and homs >= BS2_HOM_MIN
    return CriterionResult(
        "BS2", name, "strong", applies=True, met=met,
        applied_strength="strong" if met else None,
        confidence="high" if homs is not None else "low",
        evidence={"gnomad_homozygotes": homs, "cutoff": BS2_HOM_MIN},
        citation=[c for c in [a.source.get("gnomad")] if c] if homs is not None else [],
        reasoning=("gnomAD homozygote count unavailable — cannot assess healthy homozygotes"
                   if homs is None
                   else f"{homs} homozygotes in gnomAD" if met
                   else f"{homs} homozygotes (below {BS2_HOM_MIN})"),
    )


@criterion("BS3")
def bs3(v: Variant, a: Annotation) -> CriterionResult:
    # The benign mirror of PS3: same evidence class (published functional assays), same reason it
    # cannot be automated — it needs a literature judgement, not a lookup.
    return _judgment(
        "BS3", "Well-established functional studies show NO damaging effect", "strong",
        "Requires literature review of functional assays — left for expert/model adjudication")


@criterion("BS4")
def bs4(v: Variant, a: Annotation) -> CriterionResult:
    return _na("BS4", "Lack of segregation in affected members of a family", "strong",
               "Requires genotypes for affected relatives — not available from a single proband VCF")


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


@criterion("BP2")
def bp2(v: Variant, a: Annotation) -> CriterionResult:
    return _na("BP2",
               "Observed in trans with a pathogenic variant (dominant gene), or in cis with one",
               "supporting",
               "Requires phasing / parental data to establish trans or cis — not available from a "
               "single proband VCF")


@criterion("BP3")
def bp3(v: Variant, a: Annotation) -> CriterionResult:
    # Needs a repeat/low-complexity track AND a "no known function" domain call. The engine carries
    # neither, and inferring "repetitive" from the sequence alone would be a guess dressed as a fact,
    # so the in-frame consequence is surfaced and the decision left explicit.
    return _judgment(
        "BP3", "In-frame indel in a repetitive region without a known function", "supporting",
        "Requires a repeat/domain annotation the engine does not carry — model adjudication",
        evidence={"consequence": v.consequence, "hgvs_p": v.hgvs_p})


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
        # PP3's identical REVEL/CADD branch cites this; BP4's did not, so a benign call could fire
        # with the Source column blank while a predictor demonstrably drove it.
        citation=[c for c in [a.source.get("insilico")] if c],
        reasoning=(f"{_insilico_names(a)} below benign cutoffs"
                   if met else ("in-silico predictors conflict — neither PP3 nor BP4 applied"
                                if direction == "conflicting"
                                else "in-silico predictors not benign / unavailable")),
    )


@criterion("BP5")
def bp5(v: Variant, a: Annotation) -> CriterionResult:
    # Whether ANOTHER finding already explains the phenotype is a report-level, cross-variant
    # judgement (and a clinical one — an alternate basis does not always exclude a second hit).
    # A per-variant evaluator cannot see the rest of the case, so this stays explicit rather than
    # being silently inferred from the candidate list.
    return _judgment(
        "BP5", "Variant found in a case with an alternate molecular basis for disease", "supporting",
        "Requires whole-case context (another finding explaining the phenotype) plus clinical "
        "judgement — left for expert/model adjudication",
        evidence={"gene": v.gene})


@criterion("BP6")
def bp6(v: Variant, a: Annotation) -> CriterionResult:
    name = "Reputable source (ClinVar) classifies the variant as benign"
    # Mirror of PP5, sharing both halves of its gate: a criteria-based (>=1-star) ClinVar Benign /
    # Likely benign assertion. "conflicting" is deliberately excluded (it starts with neither
    # token) — a conflicted record is not a reputable benign call. The variant's OWN assertion,
    # not a residue cross-match.
    met = _clinvar_says(a, "benign", "likely benign") and _clinvar_reviewed(a)
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
