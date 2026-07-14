"""Annotation detection + the ACMG capability map (Stages 3 & 5 of the guided flow)."""
from vcf2report import status
from vcf2report.inspect import _detect_annotation, analysis_capabilities
from vcf2report.models import Variant


def _v(**kw):
    return Variant(chrom="1", pos=1, ref="A", alt="T", **kw)


def test_detect_consequence_field_wins():
    a, s = _detect_annotation([_v(consequence="missense_variant")])
    assert a and s == "consequence"


def test_detect_snpeff_ann_info():
    a, s = _detect_annotation([_v(info={"ANN": "T|missense_variant|MODERATE|GENE"})])
    assert a and s == "SnpEff ANN"


def test_detect_vep_csq_info():
    a, s = _detect_annotation([_v(info={"CSQ": "T|missense_variant|MODERATE|GENE"})])
    assert a and s == "VEP CSQ"


def test_detect_unannotated():
    a, s = _detect_annotation([_v()])
    assert not a and s is None


def _insp(annotated=True, source="consequence", build="GRCh38"):
    return {"annotated": annotated, "annotation_source": source, "build": build,
            "total_variants": 100, "pass_filter": 90}


def _rd(gnomad=True, am=True, clinvar=True, hpo=True):
    def s(p):
        return {"present": p, "path": None, "enables": ""}
    return {"stores": {"gnomad_parquet": s(gnomad), "alphamissense": s(am),
                       "clinvar_tabix": s(clinvar), "hpo": s(hpo)},
            "bundled_local_data": {"clinvar_slice": clinvar}}


def test_caps_unannotated_limits_lof_criteria():
    caps = analysis_capabilities("x", inspection=_insp(annotated=False, source=None), rd=_rd())["criteria"]
    assert caps["PVS1 (LoF)"]["status"] == "limited"
    assert caps["PM4 (in-frame / stop-loss)"]["status"] == "limited"


def test_caps_single_proband_segregation_na():
    caps = analysis_capabilities("x", inspection=_insp(), rd=_rd())["criteria"]
    assert caps["PS2 / PM3 / PM6 / PP1 / BS4 (segregation)"]["status"] == "na"


def test_caps_gnomad_absent_limits_frequency():
    # The never-fabricate-absence invariant surfaced honestly: no store -> frequency limited.
    caps = analysis_capabilities("x", inspection=_insp(), rd=_rd(gnomad=False))["criteria"]
    assert caps["PM2 / BA1 / BS1 (frequency)"]["status"] == "limited"
    assert "over-call" in caps["PM2 / BA1 / BS1 (frequency)"]["reason"]


def test_caps_missense_limited_without_alphamissense():
    caps = analysis_capabilities("x", inspection=_insp(), rd=_rd(am=False))["criteria"]
    assert caps["PP3 / BP4 (missense)"]["status"] == "limited"


def test_caps_all_available_when_fully_provisioned():
    caps = analysis_capabilities("x", hpo_given=True, inspection=_insp(), rd=_rd())["criteria"]
    assert caps["PVS1 (LoF)"]["status"] == "available"
    assert caps["PP3 / BP4 (missense)"]["status"] == "available"
    assert caps["PM2 / BA1 / BS1 (frequency)"]["status"] == "available"
    assert caps["PP4 (phenotype)"]["status"] == "available"


def test_caps_pp4_na_without_hpo_store():
    caps = analysis_capabilities("x", hpo_given=True, inspection=_insp(), rd=_rd(hpo=False))["criteria"]
    assert caps["PP4 (phenotype)"]["status"] == "na"


def test_readiness_stores_carry_enables_notes():
    r = status.readiness()
    for k in ("gnomad_parquet", "alphamissense", "clinvar_tabix", "hpo"):
        assert r["stores"][k]["enables"], f"{k} must explain what it enables"
    assert r["python"] and r["package_importable"] is True
