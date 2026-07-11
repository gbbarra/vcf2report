"""Legacy SnpEff EFF parsing + DP-from-AD depth fallback (real-world exome shapes)."""
from vcf2report.vcf import annparse
from vcf2report.vcf.parse import _sample_metrics


# ---------------------------------------------------------------------------
# Legacy SnpEff EFF (pre-ANN format)
# ---------------------------------------------------------------------------
def test_eff_missense_fields():
    eff = ("DOWNSTREAM(MODIFIER||3160||178|SAMD11|protein_coding|CODING|ENST00000420190||1),"
           "NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Tgg/Cgg|W169R|540|SAMD11|protein_coding|"
           "CODING|ENST00000455979|4|1|WARNING)")
    r = annparse.parse_snpeff_eff(eff)
    assert r["gene"] == "SAMD11"
    assert r["consequence"] == "missense_variant"          # SnpEff name -> SO term
    assert r["hgvs_p"] == "p.Trp169Arg"                    # one-letter -> HGVS.p
    assert r["transcript"] == "ENST00000455979" and r["exon"] == "4"
    assert r["hgvs_c"] is None                             # EFF has no HGVS.c


def test_eff_severity_and_stop():
    # HIGH stop_gained must win over a MODERATE missense on another transcript.
    eff = ("NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|x|Q59R|100|G|protein_coding|CODING|ENST1|2|1),"
           "STOP_GAINED(HIGH|NONSENSE|x|Q59*|100|G|protein_coding|CODING|ENST2|2|1)")
    r = annparse.parse_snpeff_eff(eff)
    assert r["consequence"] == "stop_gained"
    assert r["hgvs_p"] == "p.Gln59Ter" and r["transcript"] == "ENST2"


def test_eff_multiallelic_allele_match():
    # Genotype_Number (last field) ties each effect to its ALT.
    eff = ("STOP_GAINED(HIGH|NONSENSE|x|Q1*|1|G|protein_coding|CODING|ENST_A|1|1),"
           "NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|x|A2V|1|G|protein_coding|CODING|ENST_B|1|2)")
    r0 = annparse.parse_snpeff_eff(eff, alt_index=0, n_alt=2)
    assert r0["transcript"] == "ENST_A" and r0["consequence"] == "stop_gained"
    r1 = annparse.parse_snpeff_eff(eff, alt_index=1, n_alt=2)
    assert r1["transcript"] == "ENST_B" and r1["consequence"] == "missense_variant"


def test_eff_via_extract():
    r = annparse.extract(
        {"EFF": "STOP_GAINED(HIGH|NONSENSE|x|R5*|9|BRCA1|protein_coding|CODING|ENST9|3|1)"},
        "T", n_alt=1)
    assert r["gene"] == "BRCA1" and r["consequence"] == "stop_gained"


# ---------------------------------------------------------------------------
# Depth fallback: many callers emit AD but no FORMAT/DP.
# ---------------------------------------------------------------------------
def test_depth_from_ad_when_no_dp():
    m = _sample_metrics("GT:AD:GQ", "0/1:10,23:99", 0)
    assert m["depth"] == 33                                 # sum(AD)
    assert m["allele_balance"] == round(23 / 33, 3)


def test_format_dp_preferred_over_ad():
    m = _sample_metrics("GT:AD:DP:GQ", "0/1:10,23:50:99", 0)
    assert m["depth"] == 50                                 # explicit DP wins
