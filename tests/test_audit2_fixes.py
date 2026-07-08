"""Regression tests for the new-surface audit fixes (multiallelic allele-awareness, etc.)."""
import json

from vcf2report.acmg.engine import classify
from vcf2report.annotate import annotate_variant
from vcf2report.models import Annotation, Classification, Variant
from vcf2report.phenopacket import load_phenopacket
from vcf2report.report.assemble import split_findings
from vcf2report.vcf import annparse
from vcf2report.vcf.parse import parse_vcf

# A multiallelic site, SnpEff-annotated, with a per-allele gnomAD AF array (Number=A):
# allele C is common (0.30), allele T is ultra-rare (0.00002).
MULTI_SNPEFF = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=ANN,Number=.,Type=String,Description="SnpEff">
##INFO=<ID=gnomad_AF,Number=A,Type=Float,Description="AF">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
2\t100\t.\tA\tC,T\t800\tPASS\tANN=C|missense_variant|MODERATE|GENEA|G|transcript|TX|protein_coding|1/2|c.1A>C|p.Xaa1Cys|||||,T|stop_gained|HIGH|GENEB|G|transcript|TX|protein_coding|1/2|c.1A>T|p.Xaa1Ter|||||;gnomad_AF=0.30,0.00002\tGT:DP:GQ:AD\t1/2:60:99:0,30,30
"""


def _write(tmp_path, text, name="m.vcf"):
    p = tmp_path / name
    p.write_text(text)
    return p


def test_multiallelic_gnomad_af_is_allele_aware(tmp_path):
    variants, _, _ = parse_vcf(_write(tmp_path, MULTI_SNPEFF))
    by_key = {v.key: v for v in variants}
    c, t = by_key["2-100-A-C"], by_key["2-100-A-T"]
    assert c.alt_index == 0 and t.alt_index == 1
    # SnpEff consequence is per-allele
    assert c.consequence == "missense_variant" and t.consequence == "stop_gained"
    # gnomAD AF must follow the allele, not always allele #1
    ac = annotate_variant(c, [])
    at = annotate_variant(t, [])
    assert ac.gnomad_af == 0.30       # common allele
    assert at.gnomad_af == 0.00002    # ultra-rare allele
    # ...and that drives PM2 correctly (met only for the truly-rare allele)
    pm2_c = next(x for x in classify(c, ac).criteria if x.code == "PM2")
    pm2_t = next(x for x in classify(t, at).criteria if x.code == "PM2")
    assert pm2_c.met is False
    assert pm2_t.met is True


def test_vep_multiallelic_picks_correct_allele():
    fmt = ["Allele", "Consequence", "SYMBOL", "HGVSc", "HGVSp", "PICK"]
    csq = "C|missense_variant|GENEA|c.1A>C|p.X|1,T|stop_gained|GENEB|c.1A>T|p.Y|1"
    rc = annparse.parse_vep(csq, "C", fmt, ref="A")
    rt = annparse.parse_vep(csq, "T", fmt, ref="A")
    assert rc["gene"] == "GENEA" and rc["consequence"] == "missense_variant"
    assert rt["gene"] == "GENEB" and rt["consequence"] == "stop_gained"


def test_minimal_alt_trims_longest_shared_prefix():
    assert annparse._minimal_alt("A", "AT") == "T"
    assert annparse._minimal_alt("AA", "AAG") == "G"
    assert annparse._minimal_alt("A", "G") == "G"
    assert annparse._minimal_alt("AT", "A") == "-"


def test_split_keeps_phenotype_matched_benign_out_of_primary():
    def mkc(gene, tier, hpo):
        return Classification(
            variant=Variant(chrom="1", pos=1, ref="A", alt="G", gene=gene),
            annotation=Annotation(hpo_match_score=hpo), criteria=[], tier=tier, rule_path="")
    primary, secondary, other = split_findings([
        mkc("A", "Benign", 0.8),                 # phenotype-matched but benign
        mkc("B", "Likely Pathogenic", 0.0),      # unrelated P/LP -> secondary
        mkc("C", "Uncertain Significance (VUS)", 0.7),  # matched VUS -> primary
    ])
    assert {c.variant.gene for c in primary} == {"C"}
    assert {c.variant.gene for c in secondary} == {"B"}
    assert {c.variant.gene for c in other} == {"A"}


def test_phenopacket_skips_hgvs_only_and_escapes_info(tmp_path):
    pkt = {
        "subject": {"id": "S1"},
        "phenotypicFeatures": [{"type": {"id": "HP:0001250"}}],
        "interpretations": [{"diagnosis": {"genomicInterpretations": [
            {"variantInterpretation": {"variationDescriptor": {
                "geneContext": {"symbol": "GENE;X"},   # ';' must be escaped
                "vcfRecord": {"chrom": "2", "pos": 100, "ref": "A", "alt": "T"},
                "allelicState": {"id": "GENO:0000135"}}}},
            {"variantInterpretation": {"variationDescriptor": {
                "geneContext": {"symbol": "GENE2"},
                "expressions": [{"syntax": "hgvs.c", "value": "NM_1:c.1A>T"}]}}},  # no vcfRecord
        ]}}],
    }
    p = tmp_path / "p.json"
    p.write_text(json.dumps(pkt))
    data = load_phenopacket(p)
    assert len(data["variants"]) == 1          # HGVS-only one skipped
    assert data["skipped_variants"] == 1
    from vcf2report.phenopacket import write_inputs
    vcf, hpo = tmp_path / "o.vcf", tmp_path / "o.hpo.txt"
    write_inputs(data, vcf, hpo)
    variants, _, _ = parse_vcf(vcf)            # must still parse (INFO escaped)
    assert len(variants) == 1
    assert ";" not in variants[0].gene         # escaped
