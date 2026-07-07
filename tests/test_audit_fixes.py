"""Regression tests for the adversarial-audit fixes."""
from vcf2report import config
from vcf2report.acmg.engine import classify
from vcf2report.annotate import annotate_variant
from vcf2report.models import Annotation, Variant
from vcf2report.pipeline import run_pipeline
from vcf2report.report.render import render_markdown
from vcf2report.vcf.parse import detect_build, parse_vcf, zygosity
from vcf2report.vcf.qc import apply_qc

MULTI = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=GENE,Number=1,Type=String,Description="g">
##INFO=<ID=CSQ,Number=1,Type=String,Description="c">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
2\t100\t.\tA\tG,T\t800\tPASS\tGENE=X;CSQ=missense_variant\tGT:DP:GQ:AD\t1/2:40:99:0,20,20
3\t200\t.\tC\tG\t800\tPASS\tGENE=Y;CSQ=missense_variant\tGT:DP:GQ:AD\t0/0:40:99:40,0
4\t300\t.\tT\tA\t800\tPASS\tGENE=Z;CSQ=missense_variant\tGT:DP:GQ:AD\t./.:40:99:20,20
"""


def _write(tmp_path, text):
    p = tmp_path / "m.vcf"
    p.write_text(text)
    return p


# --- #11/#14 multiallelic per-allele metrics ---------------------------------
def test_multiallelic_per_allele_metrics(tmp_path):
    variants, build, _ = parse_vcf(_write(tmp_path, MULTI))
    assert build == "GRCh38"
    by_key = {v.key: v for v in variants}
    g, t = by_key["2-100-A-G"], by_key["2-100-A-T"]
    # 1/2 compound het -> each split allele is het, each with its own AD-based AB
    assert g.zygosity == "het" and t.zygosity == "het"
    assert g.allele_balance == 0.5 and t.allele_balance == 0.5
    # non-carriers keep zygosity None
    assert by_key["3-200-C-G"].zygosity is None   # 0/0 hom-ref
    assert by_key["4-300-T-A"].zygosity is None   # ./. no-call


# --- #7 carrier gate ---------------------------------------------------------
def test_carrier_gate_drops_noncarriers(tmp_path):
    variants, _, _ = parse_vcf(_write(tmp_path, MULTI))
    kept, dropped = apply_qc(variants)
    dropped_keys = {v.key for v, _ in dropped}
    assert "3-200-C-G" in dropped_keys   # hom-ref never reportable
    assert "4-300-T-A" in dropped_keys   # no-call never reportable
    assert "2-100-A-G" in {v.key for v in kept}


# --- #1/#4 build detection ---------------------------------------------------
def test_detect_build_does_not_misread_grch37_as_grch38():
    hdr = ["##reference=/ref/human_g1k_v37.fasta",
           "##contig=<ID=7,length=159138663>"]   # length contains "38"
    assert detect_build(hdr) == "GRCh37"


def test_detect_build_unknown_returns_none():
    assert detect_build(["##contig=<ID=1,length=248956422>"]) is None


# --- #2/#5 shared zygosity helper -------------------------------------------
def test_zygosity_helper_edge_cases():
    assert zygosity(["1", "2"], 1) == "het"    # compound het, allele 1
    assert zygosity(["1", "2"], 2) == "het"    # compound het, allele 2
    assert zygosity([".", "."], 1) is None     # no-call
    assert zygosity(["-1", "-1"], 1) is None   # cyvcf2 missing encoding
    assert zygosity(["0", "0"], 1) is None      # hom-ref
    assert zygosity(["1", "1"], 1) == "hom"
    assert zygosity(["0", "1"], 1) == "het"
    assert zygosity(["0", "2"], 1) is None      # carries a different ALT


# --- #8 build-gated annotation ----------------------------------------------
def test_untrusted_build_skips_coordinate_annotation():
    v = Variant(chrom="2", pos=166003360, ref="C", alt="T", gene="SCN1A",
                consequence="stop_gained")
    a = annotate_variant(v, [], build_trusted=False)
    assert a.gnomad_af is None            # unknown, not a fabricated 0.0
    assert a.clinvar_significance is None
    # PM2 must not fire on unavailable frequency
    pm2 = next(c for c in classify(v, a).criteria if c.code == "PM2")
    assert pm2.met is False


# --- #12 PP3/BP4 mutual exclusivity -----------------------------------------
def test_insilico_conflict_fires_neither_pp3_nor_bp4():
    v = Variant(chrom="1", pos=1, ref="A", alt="G", gene="X",
                consequence="missense_variant")
    a = Annotation(gnomad_af=0.0, abraom_af=0.0, revel=0.95, cadd_phred=5.0)
    codes = {c.code: c for c in classify(v, a).criteria}
    assert codes["PP3"].met is False      # revel says pathogenic...
    assert codes["BP4"].met is False      # ...cadd says benign -> conflict -> neither


# --- #3/#10 report render ----------------------------------------------------
def test_report_render_no_placeholder_and_expected_content():
    hpo = ["HP:0001250", "HP:0001263", "HP:0002133"]
    report = run_pipeline(config.SAMPLE_VCF, hpo_terms=hpo)
    md = render_markdown(report)
    assert "(timestamp filled by caller)" not in md
    assert report.generated                       # real ISO timestamp
    for gene in ("SCN1A", "KCNQ2", "CACNA1A"):
        assert gene in md
    assert "ABraOM" in md and "OBSCN" in md        # differentiator callout
    assert "Candidates classified: 3" in md
