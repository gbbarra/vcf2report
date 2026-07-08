"""Validate the cyvcf2 (htslib) parse path and that it agrees with the pure reader.

cyvcf2 is optional; skipped if unavailable. Uses a proper-header VCF (cyvcf2
requires declared contigs/FORMATs), covering multiallelic, star alleles, and
missing DP/GQ sentinels.
"""
import pytest

pytest.importorskip("cyvcf2")

from vcf2report.vcf.parse import parse_vcf  # noqa: E402

PROPER = """##fileformat=VCFv4.2
##reference=GRCh38
##contig=<ID=2,length=242193529>
##INFO=<ID=GENE,Number=1,Type=String,Description="g">
##INFO=<ID=CSQ,Number=1,Type=String,Description="c">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS
2\t166003360\t.\tC\tT\t800\tPASS\tGENE=SCN1A;CSQ=stop_gained\tGT:DP:GQ:AD\t0/1:45:99:22,23
2\t100\t.\tA\tG,T\t800\tPASS\tGENE=X;CSQ=missense_variant\tGT:DP:GQ:AD\t1/2:40:99:0,20,20
2\t200\t.\tC\tT,*\t800\tPASS\tGENE=Y;CSQ=missense_variant\tGT:DP:GQ:AD\t0/1:50:99:25,25,0
2\t300\t.\tA\tG\t800\tPASS\tGENE=Z;CSQ=missense_variant\tGT:DP:GQ:AD\t0/1:.:.:.,.
"""


def _tuples(variants):
    return sorted((v.key, v.gene, v.consequence, v.zygosity, v.depth, v.gq,
                   v.allele_balance) for v in variants)


def test_cyvcf2_path_and_agreement_with_pure(tmp_path, monkeypatch):
    p = tmp_path / "p.vcf"
    p.write_text(PROPER)

    monkeypatch.delenv("VCF2REPORT_NO_CYVCF2", raising=False)   # enable cyvcf2
    cy = parse_vcf(p)[0]
    monkeypatch.setenv("VCF2REPORT_NO_CYVCF2", "1")             # force pure
    pure = parse_vcf(p)[0]

    keys = {v.key for v in cy}
    assert "2-200-C-*" not in keys              # star allele skipped
    assert {"2-166003360-C-T", "2-100-A-G", "2-100-A-T", "2-200-C-T", "2-300-A-G"} <= keys

    scn = next(v for v in cy if v.key == "2-166003360-C-T")
    assert scn.depth == 45 and scn.gq == 99 and scn.zygosity == "het"

    v300 = next(v for v in cy if v.key == "2-300-A-G")          # DP/GQ '.' sentinels
    assert v300.depth is None and v300.gq is None

    # cyvcf2 and the pure reader must agree on the essentials.
    assert _tuples(cy) == _tuples(pure)
