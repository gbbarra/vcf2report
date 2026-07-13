"""gnomAD DuckDB/Parquet source: batch prime, the full/partial absence guard, case."""
import json

import pytest

from vcf2report import config
from vcf2report.annotate import gnomad, gnomad_parquet
from vcf2report.models import Variant

duckdb = pytest.importorskip("duckdb")


def _make_parquet(tmp_path, mode=None, contigs=("chr1",)):
    p = tmp_path / "g.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 0.30, 0.50, 300, 1000, 20, 0.48, 'nfe'),
        ('chr1', 200, 'C', 'G', 0.001, 0.002, 1, 1000, 0, 0.0005, 'afr'))
        AS t(chrom, pos, ref, alt, af, af_grpmax, ac, an, nhomalt, faf95, grpmax_pop))
        TO '{p}' (FORMAT PARQUET)""")
    con.close()
    if mode:   # a build that vouches for the store writes this sidecar
        (tmp_path / "g.parquet.meta.json").write_text(
            json.dumps({"mode": mode, "contigs": list(contigs)}))
    return p


@pytest.fixture
def parquet(tmp_path, monkeypatch):
    def _setup(mode=None, contigs=("chr1",)):
        p = _make_parquet(tmp_path, mode, contigs)
        monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
        gnomad_parquet._reset_for_tests()
        return p
    yield _setup
    gnomad_parquet._reset_for_tests()


def _v(chrom, pos, ref, alt):
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


def test_prime_and_match(parquet):
    parquet()
    n = gnomad_parquet.prime([_v("1", 100, "A", "T"), _v("chr1", 200, "C", "G")])
    assert n == 2
    r = gnomad_parquet.get("1-100-A-T")
    # af is the grpmax (popmax) value the ACMG engine cites, not the overall af.
    assert r["af"] == 0.50 and r["faf95"] == 0.48 and r["pop"] == "nfe" and r["hom"] == 20


def test_partial_default_never_fabricates_absence(parquet):
    # No sidecar -> partial -> a covered-site miss must be None (fall through), NOT 0.0.
    parquet()
    gnomad_parquet.prime([_v("1", 999, "A", "T")])
    assert gnomad_parquet.get("1-999-A-T") is None


def test_full_covered_absence_is_zero(parquet):
    # A build that declares mode=full + contigs may assert absence on those contigs.
    parquet(mode="full", contigs=("chr1",))
    gnomad_parquet.prime([_v("1", 999, "A", "T")])
    assert gnomad_parquet.get("1-999-A-T") == {
        "af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}


def test_bed_mode_asserts_absence_only_inside_intervals(tmp_path, monkeypatch):
    # A panel (mode='bed') store asserts absence ONLY inside a covered BED interval —
    # sound because the store is complete there; off-panel it stays unprimed (no fake 0.0).
    p = _make_parquet(tmp_path)
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t50\t250\n")     # 0-based half-open -> covers 1-based 51..250
    (tmp_path / "g.parquet.meta.json").write_text(json.dumps(
        {"mode": "bed", "contigs": ["chr1"], "bed_path": str(bed)}))
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
    gnomad_parquet._reset_for_tests()
    gnomad_parquet.prime([_v("1", 100, "A", "T"),    # matches parquet -> served
                          _v("1", 150, "A", "T"),    # in BED, absent -> genuine 0.0
                          _v("1", 5000, "A", "T")])   # off BED, absent -> unprimed
    assert gnomad_parquet.get("1-100-A-T")["af"] == 0.50
    assert gnomad_parquet.get("1-150-A-T") == {
        "af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
    assert gnomad_parquet.get("1-5000-A-T") is None
    gnomad_parquet._reset_for_tests()


def test_bed_mode_without_bed_file_stays_partial(tmp_path, monkeypatch):
    # mode='bed' but the BED can't be loaded -> fall back to partial (never assert absence).
    p = _make_parquet(tmp_path)
    (tmp_path / "g.parquet.meta.json").write_text(json.dumps(
        {"mode": "bed", "contigs": ["chr1"], "bed_path": str(tmp_path / "missing.bed")}))
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
    gnomad_parquet._reset_for_tests()
    gnomad_parquet.prime([_v("1", 150, "A", "T")])
    assert gnomad_parquet.get("1-150-A-T") is None      # no BED -> safe, unprimed
    gnomad_parquet._reset_for_tests()


def test_bed_mode_non_pass_variant_is_not_a_false_absence(tmp_path, monkeypatch):
    # Regression: a non-PASS record INSIDE the covered BED must not be turned into a fake
    # absence (which would fire a spurious PM2 -> over-call, e.g. a 99.98% AS_VQSR variant on
    # a healthy exome). It is present -> served af=None; a truly-absent in-BED locus -> af=0.
    p = tmp_path / "g.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 'PASS',    0.30, 0.50, 300, 1000, 20, 0.48, 'nfe'),
        ('chr1', 120, 'G', 'GA', 'AS_VQSR', 0.9998, 0.9998, 999, 1000, 400, 0.99, 'nfe'))
        AS t(chrom, pos, ref, alt, filter, af, af_grpmax, ac, an, nhomalt, faf95, grpmax_pop))
        TO '{p}' (FORMAT PARQUET)""")
    con.close()
    bed = tmp_path / "panel.bed"
    bed.write_text("chr1\t50\t250\n")
    (tmp_path / "g.parquet.meta.json").write_text(json.dumps(
        {"mode": "bed", "contigs": ["chr1"], "bed_path": str(bed)}))
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
    gnomad_parquet._reset_for_tests()
    gnomad_parquet.prime([_v("1", 120, "G", "GA"), _v("1", 130, "A", "T")])
    rec = gnomad_parquet.get("1-120-G-GA")                 # present but non-PASS
    assert rec is not None and rec["af"] is None           # NOT a fake absence (af != 0.0)
    assert gnomad_parquet.get("1-130-A-T") == {            # truly absent in-BED -> real absence
        "af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
    gnomad_parquet._reset_for_tests()


def test_full_uncovered_contig_left_unprimed(parquet):
    # chr2 isn't in the full store's contigs -> unprimed -> caller falls back.
    parquet(mode="full", contigs=("chr1",))
    gnomad_parquet.prime([_v("2", 500, "A", "T")])
    assert gnomad_parquet.get("2-500-A-T") is None


def test_case_insensitive_alleles(parquet):
    # Lowercase input alleles must still match (VCF alleles are case-insensitive) —
    # a byte-exact join would fabricate an absence. Variant.key upper-cases, so the
    # lowercase input resolves to the uppercase gnomAD row.
    parquet(mode="full", contigs=("chr1",))
    v = _v("1", 100, "a", "t")
    gnomad_parquet.prime([v])
    assert gnomad_parquet.get(v.key)["af"] == 0.50


def test_non_pass_variant_present_but_af_unavailable(tmp_path, monkeypatch):
    # A non-PASS gnomAD record (AS_VQSR/InbreedingCoeff/AC0 artifact) must NOT be served as an
    # authoritative frequency — its AF/faf95 would fire BA1/BS1 and mask a real pathogenic
    # variant (ClinGen/Whiffin filtering-AF is PASS-only). But the variant IS present, so it must
    # NOT be turned into a fake absence either (that would fire a spurious PM2 -> over-call).
    # Correct behaviour: served with af=None (present, frequency unavailable).
    p = tmp_path / "g.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 'PASS',    0.30, 0.50, 300, 1000, 20, 0.48, 'nfe'),
        ('chr1', 200, 'C', 'G', 'AS_VQSR', 0.40, 0.45, 400, 1000,  0, 0.44, 'afr'))
        AS t(chrom, pos, ref, alt, filter, af, af_grpmax, ac, an, nhomalt, faf95, grpmax_pop))
        TO '{p}' (FORMAT PARQUET)""")
    con.close()
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
    gnomad_parquet._reset_for_tests()
    gnomad_parquet.prime([_v("1", 100, "A", "T"), _v("1", 200, "C", "G")])
    assert gnomad_parquet.get("1-100-A-T")["af"] == 0.50   # PASS -> served
    rec = gnomad_parquet.get("1-200-C-G")                   # non-PASS artifact
    assert rec is not None and rec["af"] is None            # present, but AF not served (not absent)
    gnomad_parquet._reset_for_tests()


def test_chrom_prefix_tolerated(parquet):
    parquet()
    gnomad_parquet.prime([_v("chr1", 100, "A", "T")])
    assert gnomad_parquet.get("1-100-A-T")["af"] == 0.50


def test_off_when_unconfigured(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "GNOMAD_PARQUET", None)
    gnomad_parquet._reset_for_tests()
    assert gnomad_parquet.available() is False
    assert gnomad_parquet.prime([_v("1", 100, "A", "T")]) == 0
    gnomad_parquet._reset_for_tests()


def test_auto_detects_local_store(tmp_path, monkeypatch):
    # With no env var, the LOCAL default store is picked up automatically once it exists
    # (so a built/fetched data/gnomad/gnomad_parquet/ is "always used" without config).
    monkeypatch.delenv("VCF2REPORT_GNOMAD_PARQUET", raising=False)
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DEFAULT_GNOMAD_PARQUET", tmp_path / "gnomad" / "gnomad_parquet")
    assert config._resolve_gnomad_parquet() is None      # not built yet -> feature off
    d = tmp_path / "gnomad" / "gnomad_parquet"
    d.mkdir(parents=True)
    assert config._resolve_gnomad_parquet() == str(d)    # present -> auto-detected


def test_env_var_overrides_local_store(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DEFAULT_GNOMAD_PARQUET", tmp_path / "gnomad" / "gnomad_parquet")
    (tmp_path / "gnomad" / "gnomad_parquet").mkdir(parents=True)
    monkeypatch.setenv("VCF2REPORT_GNOMAD_PARQUET", "/somewhere/else.parquet")
    assert config._resolve_gnomad_parquet() == "/somewhere/else.parquet"


def test_lookup_prefers_primed_parquet(parquet, monkeypatch):
    from vcf2report.annotate import cache
    parquet()
    monkeypatch.setattr(cache, "get", lambda *a, **k: None)
    gnomad_parquet.prime([_v("1", 100, "A", "T")])
    r = gnomad.lookup(_v("1", 100, "A", "T"))
    assert r["af"] == 0.50 and "parquet" in r["_source"]


def test_parquet_wins_over_stale_cache(parquet, monkeypatch):
    # M1: a stale/wrong disk-cache entry must NOT shadow the fresh parquet answer —
    # the fresh authoritative source is checked before the persisted cache.
    from vcf2report.annotate import cache
    parquet()
    gnomad_parquet.prime([_v("1", 100, "A", "T")])
    monkeypatch.setattr(cache, "get",
                        lambda *a, **k: {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None})
    r = gnomad.lookup(_v("1", 100, "A", "T"))
    assert r["af"] == 0.50 and "parquet" in r["_source"]
