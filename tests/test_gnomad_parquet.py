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


def test_prefers_generic_build_over_imported(tmp_path, monkeypatch):
    # The vcf2report exome build wins over an imported store; neither is deleted.
    monkeypatch.delenv("VCF2REPORT_GNOMAD_PARQUET", raising=False)
    monkeypatch.setattr(config, "DATA_DIR", tmp_path)
    monkeypatch.setattr(config, "DEFAULT_GNOMAD_PARQUET", tmp_path / "gnomad" / "gnomad_parquet")
    (tmp_path / "gnomad" / "gnomad_parquet").mkdir(parents=True)
    (tmp_path / "gnomad" / "gnomad_parquet_generic").mkdir(parents=True)
    assert config._resolve_gnomad_parquet() == str(tmp_path / "gnomad" / "gnomad_parquet_generic")


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
