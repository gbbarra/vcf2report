"""VCF parsing + QC tests against the bundled sample."""
from vcf2report import config
from vcf2report.vcf.parse import parse_vcf, detect_build
from vcf2report.vcf.qc import apply_qc


def test_sample_parses():
    variants, build, header = parse_vcf(config.SAMPLE_VCF)
    assert build == "GRCh38"
    assert len(variants) == 10
    scn1a = next(v for v in variants if v.gene == "SCN1A")
    assert scn1a.key == "2-166003360-C-T"
    assert scn1a.consequence == "stop_gained"
    assert scn1a.is_lof is True
    assert scn1a.zygosity == "het"


def test_qc_drops_low_depth_and_lowqual():
    variants, _, _ = parse_vcf(config.SAMPLE_VCF)
    kept, dropped = apply_qc(variants)
    dropped_genes = {v.gene for v, _ in dropped}
    assert "CFTR" in dropped_genes   # DP=6 < 10
    assert "HBB" in dropped_genes    # FILTER=LowQual
    assert len(kept) == 8


def test_variant_key_strips_chr_prefix():
    from vcf2report.models import Variant
    v = Variant(chrom="chr2", pos=100, ref="A", alt="T")
    assert v.key == "2-100-A-T"


def test_detect_build():
    assert detect_build(["##reference=GRCh37"]) == "GRCh37"
