"""Golden ACMG classification tests — the differentiator's safety net."""
from vcf2report.acmg.engine import classify
from vcf2report.acmg import rules
from vcf2report.models import Annotation, CriterionResult, Variant


def _lof_variant():
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A",
                consequence="stop_gained", hgvs_p="p.Arg612Ter", zygosity="het")
    a = Annotation(gene_lof_intolerant=True, gnomad_af=0.0, abraom_af=0.0,
                   clinvar_significance="Pathogenic", clinvar_accession="VCV000012345",
                   hpo_match_score=1.0, hpo_matched_terms=["HP:0001250"])
    return v, a


def test_lof_absent_phenotype_match_is_pathogenic():
    tier = classify(*_lof_variant()).tier
    assert tier == "Pathogenic"


def test_common_variant_is_benign():
    v = Variant(chrom="1", pos=1, ref="A", alt="G", gene="X", consequence="missense_variant")
    a = Annotation(gnomad_af=0.12, abraom_af=0.10)
    assert classify(v, a).tier == "Benign"


def test_rare_missense_no_other_evidence_is_vus():
    v = Variant(chrom="3", pos=3, ref="A", alt="G", gene="Y", consequence="missense_variant")
    a = Annotation(gnomad_af=0.0, abraom_af=0.0)
    assert "VUS" in classify(v, a).tier


def test_abraom_blocks_pm2():
    """A variant absent from gnomAD but present in ABraOM must NOT earn PM2."""
    v = Variant(chrom="1", pos=228208000, ref="G", alt="A", gene="OBSCN",
                consequence="missense_variant")
    a = Annotation(gnomad_af=0.0, abraom_af=0.03)
    result = classify(v, a)
    pm2 = next(c for c in result.criteria if c.code == "PM2")
    assert pm2.met is False, "ABraOM presence should block PM2"


def test_pm2_not_met_when_gnomad_frequency_unknown():
    """Unknown gnomAD AF (None) must not be treated as absence -> PM2 not met."""
    v = Variant(chrom="1", pos=228208000, ref="G", alt="A", gene="OBSCN",
                consequence="missense_variant")
    a = Annotation(gnomad_af=None, abraom_af=0.0)  # frequency unavailable
    result = classify(v, a)
    pm2 = next(c for c in result.criteria if c.code == "PM2")
    assert pm2.met is False
    assert pm2.confidence == "low"


def test_conflicting_evidence_is_vus():
    """Pathogenic + benign evidence resolves to VUS, not a forced call."""
    crits = [
        CriterionResult("PVS1", "n", "very_strong", applies=True, met=True,
                        applied_strength="very_strong"),
        CriterionResult("PM2", "n", "moderate", applies=True, met=True,
                        applied_strength="moderate"),
        CriterionResult("BA1", "n", "stand_alone", applies=True, met=True,
                        applied_strength="stand_alone"),
    ]
    tier, path = rules.combine(crits)
    assert "VUS" in tier
    assert "conflicting" in path


def test_combining_rule_lp_from_one_strong_one_moderate():
    crits = [
        CriterionResult("PS1", "n", "strong", applies=True, met=True, applied_strength="strong"),
        CriterionResult("PM2", "n", "moderate", applies=True, met=True, applied_strength="moderate"),
    ]
    tier, _ = rules.combine(crits)
    assert tier == "Likely Pathogenic"
