"""Recessive carrier alleles must not be presented as diagnostic findings.

Opening PVS1 to recessive genes (the constraint gate is blind to them — carriers are healthy,
so the gene never looks constrained) makes the engine call the het null alleles that every
healthy person carries. The Pathogenic TIER is correct: ACMG classifies the variant, not the
patient. What is NOT correct is letting that variant-level tier become a patient-level claim.

Two failures this guards, both measured on the real annotated cohort before the fix:
  * a het LIPA/SKIC2 carrier reached PRIMARY and DISPLACED the planted COPZ1 into "other",
    with the conclusion asserting "Likely explanatory finding: LIPA — Pathogenic";
  * a het PMM2 carrier was reported co-equal with the true POGZ diagnosis.
Phenotype routing cannot save this: recessive disease genes have exactly the phenotypes a
proband presents with, so the carrier clears HPO_RELATED_MIN too.
"""
from vcf2report.models import Annotation, Classification, Variant
from vcf2report.report.assemble import (carrier_findings, is_unconfirmed_ar_carrier,
                                        split_findings)
import pytest

from vcf2report.annotate import inheritance


@pytest.fixture
def moi(monkeypatch):
    def _set(mapping):
        monkeypatch.setattr(inheritance, "_gene_moi",
                            {g.upper(): frozenset(v) for g, v in mapping.items()})
    return _set


def _c(gene, tier="Pathogenic", zyg="het", hpo=0.9, pos=1):
    return Classification(
        variant=Variant(chrom="1", pos=pos, ref="A", alt="G", gene=gene, zygosity=zyg),
        annotation=Annotation(hpo_match_score=hpo, hpo_best_match=hpo),
        criteria=[], tier=tier, rule_path="")


def test_lone_het_in_recessive_gene_never_reaches_primary(moi):
    """Even when it matches the phenotype perfectly — which it will, since the gene's
    disease IS the kind of phenotype the proband has."""
    moi({"RECGENE": ["AR"], "DOMGENE": ["AD"]})
    primary, secondary, other = split_findings([_c("RECGENE"), _c("DOMGENE")])
    assert {c.variant.gene for c in primary} == {"DOMGENE"}
    assert {c.variant.gene for c in other} == {"RECGENE"}


def test_carrier_does_not_displace_the_true_diagnosis(moi):
    """The measured failure: POGZ (AD, the real answer) reported co-equal with a PMM2 carrier."""
    moi({"POGZ": ["AD"], "PMM2": ["AR"]})
    primary, _sec, other = split_findings([_c("POGZ", tier="Likely Pathogenic"), _c("PMM2")])
    assert [c.variant.gene for c in primary] == ["POGZ"]
    assert [c.variant.gene for c in other] == ["PMM2"]


def test_homozygous_in_recessive_gene_stays_diagnostic(moi):
    """Biallelic is exactly the diagnostic genotype — the fix must never hide it."""
    moi({"RECGENE": ["AR"]})
    primary, _sec, _other = split_findings([_c("RECGENE", zyg="hom")])
    assert [c.variant.gene for c in primary] == ["RECGENE"]


def test_two_hits_in_recessive_gene_stay_diagnostic(moi):
    """Possible compound heterozygote: we cannot phase it, but the clinician must see it."""
    moi({"RECGENE": ["AR"]})
    primary, _sec, _other = split_findings([_c("RECGENE", pos=1), _c("RECGENE", pos=2)])
    assert len(primary) == 2


def test_het_in_gene_with_any_dominant_disease_stays_diagnostic(moi):
    """ATM is AR (ataxia-telangiectasia) AND AD (cancer risk): a het there can be diagnostic."""
    moi({"ATM": ["AD", "AR"]})
    primary, _sec, _other = split_findings([_c("ATM")])
    assert [c.variant.gene for c in primary] == ["ATM"]


def test_recessive_sf_carrier_is_not_an_actionable_secondary(moi):
    """ACMG SF v3.2: the recessive SF genes (ATP7B, MUTYH, BTD, GAA, HFE...) are reportable
    ONLY when biallelic. Telling a healthy 1-in-90 Wilson-disease carrier they have an
    actionable ATP7B finding is a guideline violation, not a conservative extra."""
    moi({"ATP7B": ["AR"], "MUTYH": ["AR"]})
    _pri, secondary, other = split_findings([_c("ATP7B", hpo=0.0), _c("MUTYH", hpo=0.0)])
    assert secondary == []
    assert {c.variant.gene for c in other} == {"ATP7B", "MUTYH"}


def test_biallelic_sf_gene_is_still_reported(moi):
    moi({"ATP7B": ["AR"]})
    _pri, secondary, _other = split_findings([_c("ATP7B", zyg="hom", hpo=0.0)])
    assert [c.variant.gene for c in secondary] == ["ATP7B"]


def test_vus_carrier_is_not_flagged(moi):
    """Only P/LP calls are carrier findings; a VUS is not a carrier claim."""
    moi({"RECGENE": ["AR"]})
    assert not is_unconfirmed_ar_carrier(_c("RECGENE", tier="Uncertain Significance (VUS)"), [])


def test_carrier_findings_surfaces_them_for_the_report(moi):
    """Routed out of the diagnosis, NOT discarded: carrier status has reproductive relevance."""
    moi({"RECGENE": ["AR"], "DOMGENE": ["AD"]})
    cls = [_c("RECGENE"), _c("DOMGENE")]
    assert [c.variant.gene for c in carrier_findings(cls)] == ["RECGENE"]


def test_conclusion_calls_a_carrier_a_carrier(moi):
    """"Clinical relevance is uncertain" is WRONG for a carrier: the relevance is known, and
    it is reproductive rather than diagnostic. Lumping it in with genuine incidental P/LP
    invites the reader to weigh a carrier allele as a diagnostic candidate."""
    from vcf2report.models import QCSummary
    from vcf2report.report.assemble import ReportModel, summarize
    moi({"RECGENE": ["AR"], "INCGENE": ["AD"]})
    rep = ReportModel(sample_id="S", hpo_terms=[], qc=QCSummary(),
                      classifications=[_c("RECGENE", hpo=0.0), _c("INCGENE", hpo=0.0)])
    text = " ".join(summarize(rep))
    assert "Carrier finding" in text
    assert "does NOT explain the indication" in text
    # the carrier must not also appear in the "relevance uncertain" incidental bullet
    inc = [l for l in summarize(rep) if l.startswith("Additional")]
    assert inc and "INCGENE" in inc[0] and "RECGENE" not in inc[0]
