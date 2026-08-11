"""Merge of per-chromosome results.json (the chunked WGS path)."""
from vcf2report.report.merge import merge_results, merge_seq_quality


def test_seq_quality_weights_by_count_not_naive_mean():
    """A 9000-variant chunk at 40x and a 1000-variant chunk at 10x must merge to
    ~37x (weighted), never 25x (naive mean of the two means)."""
    a = {"n_variants": 9000, "n_with_dp": 9000, "dp_mean": 40.0, "dp_median": 40.0,
         "n_snv": 8000, "n_het": 5000, "n_hom": 4000, "n_indel": 1000,
         "n_sites": 9000, "n_multiallelic_sites": 90, "n_het_ab": 5000,
         "pct_het_ab_balanced": 98.0}
    b = {"n_variants": 1000, "n_with_dp": 1000, "dp_mean": 10.0, "dp_median": 10.0,
         "n_snv": 900, "n_het": 600, "n_hom": 400, "n_indel": 100,
         "n_sites": 1000, "n_multiallelic_sites": 10, "n_het_ab": 600,
         "pct_het_ab_balanced": 90.0}
    m = merge_seq_quality([a, b])
    assert m["n_variants"] == 10000
    assert m["dp_mean"] == 37.0            # (40*9000 + 10*1000)/10000, not 25.0
    # ratios recomputed exactly from summed counts
    assert m["het_hom_ratio"] == round(5600 / 4400, 2)
    assert m["indel_snv_ratio"] == round(1100 / 8900, 3)
    assert m["pct_multiallelic"] == round(100 * 100 / 10000, 4)


def test_merge_concatenates_findings_and_sums_funnel():
    r1 = {"sample_id": "S", "build": "GRCh38",
          "qc": {"total_variants": 9000, "after_qc": 8000, "candidates": 1,
                 "qc_rescued": [{"g": "A"}]},
          "seq_quality": {"n_variants": 9000, "n_with_dp": 9000, "dp_mean": 30.0},
          "classifications": [{"gene": "A"}],
          "buckets": {"primary": [{"gene": "A"}], "secondary": []},
          "clinvar_do_not_dismiss": [], "conclusion": ["c1"],
          "timings": {"total_s": 2.0}, "generated": "2026-08-11T10:00:00+00:00"}
    r2 = {"sample_id": "S", "build": "GRCh38",
          "qc": {"total_variants": 1000, "after_qc": 900, "candidates": 2,
                 "qc_rescued": []},
          "seq_quality": {"n_variants": 1000, "n_with_dp": 1000, "dp_mean": 20.0},
          "classifications": [{"gene": "B"}, {"gene": "C"}],
          "buckets": {"primary": [], "secondary": [{"gene": "B"}]},
          "clinvar_do_not_dismiss": [{"gene": "C"}], "conclusion": ["c2"],
          "timings": {"total_s": 0.5}, "generated": "2026-08-11T10:05:00+00:00"}
    m = merge_results([r1, r2])
    assert m["qc"]["total_variants"] == 10000
    assert m["qc"]["candidates"] == 3
    assert m["qc"]["qc_rescued"] == [{"g": "A"}]          # list concatenated, not summed
    assert [c["gene"] for c in m["classifications"]] == ["A", "B", "C"]
    assert [c["gene"] for c in m["buckets"]["primary"]] == ["A"]
    assert [c["gene"] for c in m["buckets"]["secondary"]] == ["B"]
    assert m["clinvar_do_not_dismiss"] == [{"gene": "C"}]
    assert m["conclusion"] == ["c1", "c2"]
    assert m["timings"]["total_s"] == 2.5
    assert m["generated"] == "2026-08-11T10:05:00+00:00"  # newest
    assert m["chunks_merged"] == 2
    assert m["seq_quality"]["dp_mean"] == 29.0            # weighted (30*9000+20*1000)/10000


def test_single_chunk_is_passthrough():
    r = {"sample_id": "S", "qc": {"candidates": 1}, "classifications": [{"gene": "A"}]}
    assert merge_results([r]) == r
