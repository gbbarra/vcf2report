"""Sequencing-quality estimate from the VCF (depth/GQ at variant sites, Ti/Tv, het:hom)."""
from vcf2report import config, pipeline
from vcf2report.models import Variant
from vcf2report.report.render import render_markdown
from vcf2report.vcf import seqqc


def _v(ref="A", alt="G", depth=30, gq=50, zyg="het", pos=1):
    return Variant(chrom="1", pos=pos, ref=ref, alt=alt, depth=depth, gq=gq, zygosity=zyg)


def test_empty_is_safe():
    q = seqqc.estimate([])
    assert q.n_variants == 0
    assert q.dp_mean is None and q.titv is None and q.het_hom_ratio is None


def test_depth_and_gq_stats():
    vs = [_v(depth=d, gq=g, pos=i) for i, (d, g) in
          enumerate([(5, 10), (15, 25), (30, 60), (40, 99)])]
    q = seqqc.estimate(vs)
    assert q.n_variants == 4 and q.n_with_dp == 4
    assert q.dp_mean == 22.5 and q.dp_median == 22.5
    assert q.dp_pct_ge10 == 75.0   # 15, 30, 40
    assert q.dp_pct_ge20 == 50.0   # 30, 40
    assert q.gq_median == 42.5
    assert q.gq_pct_ge20 == 75.0   # 25, 60, 99


def test_titv_over_snvs_only():
    # transitions A>G, C>T ; transversions A>C, G>T ; plus one indel (excluded)
    vs = [_v("A", "G", pos=1), _v("C", "T", pos=2), _v("A", "C", pos=3), _v("G", "T", pos=4),
          Variant(chrom="1", pos=5, ref="AT", alt="A")]
    q = seqqc.estimate(vs)
    assert q.n_snv == 4                 # indel excluded
    assert q.titv == 1.0               # 2 transitions / 2 transversions


def test_het_hom_ratio():
    q = seqqc.estimate([_v(zyg="het", pos=1), _v(zyg="het", pos=2), _v(zyg="hom", pos=3)])
    assert q.n_het == 2 and q.n_hom == 1 and q.het_hom_ratio == 2.0


def test_assay_guess_bands():
    assert seqqc._assay_guess(600_000) == "whole-genome-scale"
    assert seqqc._assay_guess(25_000) == "exome / large-panel-scale"
    assert seqqc._assay_guess(500) == "targeted-panel-scale"
    assert seqqc._assay_guess(50) == "small / demo VCF"


def test_missing_depth_gq_degrade_to_none():
    vs = [_v(depth=None, gq=None, pos=1), _v(depth=None, gq=None, pos=2)]
    q = seqqc.estimate(vs)
    assert q.n_with_dp == 0 and q.dp_mean is None
    assert q.n_with_gq == 0 and q.gq_median is None
    assert q.n_variants == 2  # still counted


def test_pipeline_attaches_and_renders():
    r = pipeline.run_pipeline(str(config.SAMPLE_VCF))
    assert r.seq_quality is not None and r.seq_quality.n_variants > 0
    assert r.seq_quality.dp_mean is not None
    md = render_markdown(r)
    assert "Sequencing quality" in md
    assert "Depth at called sites" in md
    assert "Indel:SNV ratio" in md and "FILTER = PASS" in md
    # honest caveat present
    assert "not genome-wide breadth" in md


def test_indel_snv_ratio():
    vs = [_v("A", "G", pos=1), _v("A", "G", pos=2),
          Variant(chrom="1", pos=3, ref="AT", alt="A"),
          Variant(chrom="1", pos=4, ref="C", alt="CGG")]
    q = seqqc.estimate(vs)
    assert q.n_snv == 2 and q.n_indel == 2 and q.indel_snv_ratio == 1.0


def test_multiallelic_by_site():
    vs = [_v(pos=10, alt="G"), _v(pos=10, alt="T"), _v(pos=20, alt="G")]
    vs[0].n_alts = vs[1].n_alts = 2   # the two records at pos 10 share a multiallelic site
    q = seqqc.estimate(vs)
    assert q.n_sites == 2 and q.n_multiallelic_sites == 1 and q.pct_multiallelic == 50.0


def test_novelty_needs_a_dbsnp_annotated_vcf():
    vs = [_v(pos=i) for i in range(10)]
    vs[0].variant_id = "rs123"                      # one stray rsID
    q = seqqc.estimate(vs)
    assert q.n_with_rsid == 1 and q.pct_novel is None   # not treated as annotated
    for v in vs:
        v.variant_id = "rs1"
    vs[0].variant_id = None                          # 9/10 annotated, 1 novel
    q2 = seqqc.estimate(vs)
    assert q2.n_with_rsid == 9 and q2.pct_novel == 10.0


def test_het_allele_balance_and_pass():
    vs = [_v(zyg="het", pos=1), _v(zyg="het", pos=2), _v(zyg="het", pos=3)]
    vs[0].allele_balance = vs[1].allele_balance = 0.5   # balanced
    vs[2].allele_balance = 0.9                          # skewed
    vs[0].filter_status = vs[1].filter_status = "PASS"
    vs[2].filter_status = "LowQual"
    q = seqqc.estimate(vs)
    assert q.n_het_ab == 3 and q.pct_het_ab_balanced == 66.7
    assert q.pct_pass == 66.7
