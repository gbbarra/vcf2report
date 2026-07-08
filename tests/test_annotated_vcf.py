"""T2: consuming a real pre-annotated VCF (SnpEff ANN / VEP CSQ + INFO)."""
from vcf2report.annotate import annotate_variant, from_vcf
from vcf2report.models import Variant
from vcf2report.vcf import annparse
from vcf2report.vcf.parse import parse_vcf

SNPEFF = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff annotations">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
2\t166003360\t.\tC\tT\t800\tPASS\tANN=T|stop_gained|HIGH|SCN1A|ENSG00000144285|transcript|ENST1|protein_coding|10/26|c.1834C>T|p.Arg612Ter|||||;gnomad_AF=0.0;CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_multiple_submitters,_no_conflicts;REVEL=0.9\tGT:DP:GQ:AD\t0/1:45:99:22,23
"""

VEP = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|HGVSc|HGVSp|CANONICAL">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
20\t63446204\t.\tG\tA\t800\tPASS\tCSQ=A|missense_variant|MODERATE|KCNQ2|ENSG1|Transcript|ENST2|protein_coding|c.637C>T|p.Arg213Trp|YES;gnomad_AF=0.00001;CLNSIG=Likely_pathogenic\tGT:DP:GQ:AD\t0/1:40:99:20,20
"""


def _write(tmp_path, text, name="a.vcf"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_snpeff_ann_parsed(tmp_path):
    variants, build, _ = parse_vcf(_write(tmp_path, SNPEFF))
    v = variants[0]
    assert v.gene == "SCN1A"
    assert v.consequence == "stop_gained"
    assert v.hgvs_p == "p.Arg612Ter"
    assert v.is_lof is True


def test_vep_csq_parsed_via_header_format(tmp_path):
    variants, _, _ = parse_vcf(_write(tmp_path, VEP))
    v = variants[0]
    assert v.gene == "KCNQ2"
    assert v.consequence == "missense_variant"
    assert v.hgvs_p == "p.Arg213Trp"


def test_info_annotations_preferred_over_lookup(tmp_path):
    variants, _, _ = parse_vcf(_write(tmp_path, SNPEFF))
    a = annotate_variant(variants[0], ["HP:0001250"])
    # gnomAD AF + ClinVar came straight from INFO (no DB lookup)
    assert a.gnomad_af == 0.0
    assert a.source["gnomad"] == "VCF INFO"
    assert a.clinvar_significance == "Pathogenic"
    assert a.source["clinvar"] == "VCF INFO"
    assert a.revel == 0.9


def test_from_vcf_extract_clinvar_underscore_normalization():
    v = Variant(chrom="1", pos=1, ref="A", alt="G",
                info={"CLNSIG": "Likely_pathogenic",
                      "CLNREVSTAT": "criteria_provided,_single_submitter"})
    out = from_vcf.extract(v)
    assert out["clinvar_significance"] == "Likely pathogenic"
    assert "criteria provided" in out["clinvar_review_status"]


def test_annparse_unit_snpeff_and_vep():
    s = annparse.parse_snpeff(
        "T|missense_variant&splice_region_variant|MODERATE|BRCA1|G|transcript|"
        "TX|protein_coding|1/2|c.1A>G|p.Met1Val|||||", "T")
    assert s["gene"] == "BRCA1"
    assert s["consequence"] == "missense_variant"  # first (most severe) term
    fmt = ["Allele", "Consequence", "SYMBOL", "HGVSc", "HGVSp"]
    vp = annparse.parse_vep("A|stop_gained|TP53|c.1A>T|p.Lys1Ter", "A", fmt)
    assert vp["gene"] == "TP53" and vp["consequence"] == "stop_gained"
