"""Absence of data must never be reported as evidence.

The engine's stated invariant: "we could not check" may never render as "we checked and it was
absent/negative". A confident wrong answer is this tool's worst failure mode, and the criteria that
violated it did so at ``confidence="high"`` with a source cited — the most convincing possible form
of a fabricated observation.

Each test below pins one criterion that used to collapse an unknown into a measured value.
"""
import pytest

from vcf2report.acmg.engine import evaluate_criteria
from vcf2report.models import Annotation, Variant


def _v(consequence="missense_variant", gene="TESTG", hgvs_p="p.Arg123Cys"):
    return Variant(chrom="1", pos=100, ref="C", alt="T", gene=gene,
                   consequence=consequence, hgvs_p=hgvs_p)


def _c(code, a, v=None):
    return next(c for c in evaluate_criteria(v or _v(), a) if c.code == code)


# --- BS2: homozygote count -------------------------------------------------
def test_bs2_does_not_invent_zero_homozygotes():
    # `homs = a.gnomad_homozygotes or 0` reported "0 homozygotes (below 2)" at confidence=high,
    # citing a source that never supplied the field — including on the build-mismatch path where
    # gnomAD was explicitly skipped.
    cr = _c("BS2", Annotation(gnomad_af=0.08, gnomad_homozygotes=None,
                              source={"gnomad": "VCF INFO"}))
    assert not cr.met
    assert "unavailable" in cr.reasoning and "cannot assess" in cr.reasoning
    assert cr.evidence["gnomad_homozygotes"] is None
    assert cr.citation == [], "cited a source for a lookup that did not happen"
    assert cr.confidence == "low"


def test_bs2_still_reports_a_real_zero():
    cr = _c("BS2", Annotation(gnomad_homozygotes=0, source={"gnomad": "gnomAD v4.1"}))
    assert not cr.met
    assert "0 homozygotes" in cr.reasoning
    assert cr.confidence == "high" and cr.citation == ["gnomAD v4.1"]


def test_bs2_fires_on_enough_homozygotes():
    assert _c("BS2", Annotation(gnomad_homozygotes=2)).met


# --- PP4: phenotype match --------------------------------------------------
def test_pp4_does_not_invent_a_zero_phenotype_match():
    # "phenotype match 0.00 below 0.6" reads as a comparison that ARGUES AGAINST the variant.
    # With no HPO terms supplied, no comparison happened at all.
    cr = _c("PP4", Annotation(hpo_match_score=None))
    assert not cr.met
    assert "no phenotype comparison" in cr.reasoning
    assert cr.citation == [] and cr.confidence == "low"


def test_hpo_match_returns_none_when_no_comparison_is_possible():
    from vcf2report.annotate import hpo
    assert hpo.match("SCN1A", [])["score"] is None        # no patient terms
    assert hpo.match(None, ["HP:0001250"])["score"] is None  # no gene
    assert hpo.match("NOT_A_GENE_XYZ", ["HP:0001250"])["score"] is None  # gene not annotated


# --- BA1 / BS1: the frequency band -----------------------------------------
def test_ba1_needs_strictly_greater_than_five_percent():
    # Richards 2015 BA1 is ">5%", as are this criterion's own name and its met-wording.
    assert not _c("BA1", Annotation(gnomad_faf95=0.05)).met, "5.00% is not >5%"
    assert _c("BA1", Annotation(gnomad_faf95=0.0501)).met


def test_bs1_does_not_claim_a_common_allele_is_below_its_cutoff():
    # An 8% allele was described as "under the 0.01 BS1 cutoff" on the line directly beneath BA1
    # reporting the same number exceeded 0.05.
    a = Annotation(gnomad_faf95=0.08)
    assert _c("BA1", a).met
    bs1 = _c("BS1", a)
    assert not bs1.met
    assert "exceeds" in bs1.reasoning and "superseded by BA1" in bs1.reasoning
    assert "under the" not in bs1.reasoning


# --- residue index: per-gene coverage --------------------------------------
def test_residue_criteria_say_not_assessed_for_a_gene_outside_the_index():
    # `available` was a global "any rows loaded" flag. The shipped frozen slice covers 15 genes, so
    # ~20,000 genes were getting "checked, no match" at high confidence for a table they were
    # never in.
    from vcf2report.annotate import annotate_variant
    v = Variant(chrom="17", pos=43093464, ref="C", alt="T", gene="BRCA1",
                consequence="missense_variant", hgvs_p="p.Arg1699Trp")
    ann = annotate_variant(v, [])
    assert ann.clinvar_residue_available is False, "BRCA1 is not in the committed slice"
    for code in ("PS1", "PM5", "PM1"):
        assert "not assessed" in _c(code, ann, v).reasoning


# --- PM1: the tolerance guard must be answerable ----------------------------
def test_pm1_does_not_fire_when_the_tolerance_metric_is_unknown():
    # `bool(None)` made "unknown" behave exactly like "proven constrained", so PM1's stand-in for
    # ACMG's "without benign variation" clause silently vanished.
    hot = {"n_residues": 5, "n_changes": 12, "window": 7, "enrichment": 4.0,
           "gene_baseline": 0.05, "available": True}
    unknown = _c("PM1", Annotation(clinvar_residue_available=True, clinvar_hotspot=hot,
                                   gene_missense_tolerant=None))
    assert not unknown.met
    assert "no gnomAD missense-constraint metric" in unknown.reasoning

    known = _c("PM1", Annotation(clinvar_residue_available=True, clinvar_hotspot=hot,
                                 gene_missense_tolerant=False))
    assert known.met


# --- citations track the decision ------------------------------------------
def test_pvs1_does_not_cite_sources_when_it_did_not_fire():
    # An unmet PVS1 printed "ClinGen Dosage Sensitivity (HI=3)" in the Source column of a row whose
    # reasoning is about consequence type.
    cr = _c("PVS1", Annotation(gene_lof_intolerant=True,
                               source={"gene_constraint": "gnomAD v4.1 constraint"}),
            _v(consequence="missense_variant"))
    assert not cr.met and cr.citation == []


def test_bp4_cites_the_predictor_that_drove_it():
    # PP3's identical REVEL/CADD branch cited this; BP4's did not.
    cr = _c("BP4", Annotation(revel=0.05, cadd_phred=3.0,
                              source={"insilico": "dbNSFP v4.4"}))
    assert cr.met and cr.citation == ["dbNSFP v4.4"]


@pytest.mark.parametrize("code,ann", [
    ("PP3", Annotation(cadd_phred=35.0, revel=None)),
    ("BP4", Annotation(revel=0.05, cadd_phred=None)),
])
def test_insilico_reasons_name_only_the_predictors_that_ran(code, ann):
    # "REVEL=None, CADD=35.0" reads as a rendering failure and hides how many predictors actually
    # contributed. REVEL is missense-only, so a missing value is the norm on LoF variants.
    cr = _c(code, ann, _v(consequence="stop_gained"))
    assert "None" not in cr.reasoning


# --- ABraOM must not be discarded by a gnomAD "absent" ----------------------
def test_abraom_frequency_survives_a_gnomad_absent_result():
    """The project's stated differentiator must not be silently dropped.

    All three gnomAD backends report an absent variant as faf95=0.0 (not None), and _benign_af
    returned faf95 the moment it existed — so a variant common in Brazilians and absent from gnomAD
    lost BA1 and BS1 while the trail asserted "= 0.0000 below 0.05". Installing the local gnomAD
    store made classification strictly worse.
    """
    a = Annotation(gnomad_af=0.0, gnomad_faf95=0.0, abraom_af=0.20)
    ba1 = _c("BA1", a)
    assert ba1.met, "20% in ABraOM must reach stand-alone benign"
    assert "ABraOM" in ba1.reasoning


def test_gnomad_faf95_still_wins_when_it_is_the_higher_value():
    # faf95 remains the preferred statistic; ABraOM only takes over when it is genuinely higher.
    a = Annotation(gnomad_faf95=0.09, abraom_af=0.01)
    cr = _c("BA1", a)
    assert cr.met and "filtering AF" in cr.reasoning


def test_benign_af_still_reports_unavailable_when_nothing_was_looked_up():
    cr = _c("BA1", Annotation())
    assert not cr.met and "cannot assess" in cr.reasoning
