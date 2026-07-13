"""In-frame indels must be recognised regardless of annotator term (VEP inframe_insertion/
deletion, SnpEff disruptive_/conservative_inframe_*, or a generic inframe_indel) — both by the
impact filter (so they are not dropped before classification) and by PM4."""
from vcf2report.acmg.criteria import pm4
from vcf2report.models import Annotation, Variant
from vcf2report.vcf.filter import filter_variants, is_impactful, is_inframe_indel

INFRAME_TERMS = ("inframe_insertion", "inframe_deletion", "disruptive_inframe_deletion",
                 "conservative_inframe_insertion", "inframe_indel")


def _v(cons):
    return Variant(chrom="1", pos=100, ref="ATCTCT", alt="A", gene="G", consequence=cons)


def test_inframe_terms_recognised_by_filter_and_pm4():
    for c in INFRAME_TERMS:
        assert is_inframe_indel(c) and is_impactful(c), c
        assert pm4(_v(c), Annotation()).met, c
    # negatives
    assert not is_inframe_indel("missense_variant")
    assert not pm4(_v("missense_variant"), Annotation()).met
    assert not is_impactful("synonymous_variant")


def test_generic_inframe_indel_survives_impact_filter():
    # a rare in-frame indel with a non-VEP consequence term must NOT be dropped pre-classification
    v = _v("inframe_indel")
    cands, _funnel = filter_variants([(v, Annotation(gnomad_af=0.0))])
    assert any(vv.consequence == "inframe_indel" for vv, _ in cands)
