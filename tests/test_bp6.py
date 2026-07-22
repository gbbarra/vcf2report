"""BP6 — reputable-source (ClinVar) benign assertion: the benign mirror of PP5.

Same gate as PP5 (criteria-based ≥1-star, Supporting strength), on the benign side. Deprecated by
the ClinGen SVI but retained as a transparent, gated line so a reviewed ClinVar *benign* call
contributes symmetrically to a reviewed pathogenic one.
"""
from vcf2report.acmg import rules
from vcf2report.acmg.criteria import all_criteria
from vcf2report.models import Annotation, CriterionResult, Variant

_bp6 = all_criteria()["BP6"]


def _v(gene="TESTG"):
    return Variant(chrom="1", pos=100, ref="A", alt="G", gene=gene,
                   consequence="missense_variant", zygosity="het")


def test_bp6_met_on_reviewed_clinvar_benign():
    cr = _bp6(_v(), Annotation(
        clinvar_significance="Benign",
        clinvar_review_status="criteria provided, multiple submitters, no conflicts",
        clinvar_accession="VCV000009999"))
    assert cr.applies and cr.met
    assert cr.applied_strength == "supporting"
    assert cr.citation == ["VCV000009999"]


def test_bp6_met_on_likely_benign_underscored_status():
    # VCF-INFO path delivers underscore-delimited review status — must normalize like PP5 does.
    cr = _bp6(_v(), Annotation(clinvar_significance="Likely benign",
                               clinvar_review_status="criteria_provided,_single_submitter"))
    assert cr.met


def test_bp6_not_met_on_zero_star():
    # "no assertion criteria provided" contains the substring "criteria provided" but is 0-star.
    cr = _bp6(_v(), Annotation(clinvar_significance="Benign",
                               clinvar_review_status="no assertion criteria provided"))
    assert cr.applies and not cr.met


def test_bp6_not_met_on_pathogenic():
    cr = _bp6(_v(), Annotation(clinvar_significance="Pathogenic",
                               clinvar_review_status="criteria provided, single submitter"))
    assert not cr.met


def test_bp6_not_met_when_no_clinvar():
    assert not _bp6(_v(), Annotation()).met


def _met(code, strength="supporting"):
    return CriterionResult(code, code, strength, applies=True, met=True, applied_strength=strength)


def test_bp6_counts_as_supporting_benign_in_combine():
    # Two Supporting benign lines (BP6 + BP4) → Likely Benign (Richards LB-2). This proves BP6 is
    # routed to the benign side of the combiner, not the pathogenic side.
    tier, path = rules.combine([_met("BP6"), _met("BP4")])
    assert tier == "Likely Benign"
    assert "BP6" in path
