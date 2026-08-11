"""VCF parsing + QC tests against the bundled sample."""
from vcf2report import config
from vcf2report.vcf.parse import parse_vcf, detect_build
from vcf2report.vcf.qc import apply_qc


def test_qual_column_is_captured(tmp_path):
    """QUAL (col 6) flows onto the Variant — '.' becomes None, a number becomes float."""
    vcf = tmp_path / "q.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
        "1\t100\t.\tA\tG\t742.5\tPASS\t.\tGT:DP:GQ\t0/1:30:99\n"
        "1\t200\t.\tC\tT\t.\tPASS\t.\tGT:DP:GQ\t0/1:30:99\n"
    )
    variants, _, _ = parse_vcf(vcf)
    by_pos = {v.pos: v for v in variants}
    assert by_pos[100].qual == 742.5
    assert by_pos[200].qual is None            # '.' → None, not 0.0
    assert by_pos[100].to_dict()["qual"] == 742.5  # reaches results.json


def test_sample_parses():
    variants, build, header = parse_vcf(config.SAMPLE_VCF)
    assert build == "GRCh38"
    assert len(variants) == 11
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
    assert len(kept) == 9


def test_variant_key_strips_chr_prefix():
    from vcf2report.models import Variant
    v = Variant(chrom="chr2", pos=100, ref="A", alt="T")
    assert v.key == "2-100-A-T"


def test_detect_build():
    assert detect_build(["##reference=GRCh37"]) == "GRCh37"
