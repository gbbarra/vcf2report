"""gnomAD DuckDB/Parquet source: batch prime + the covered-contig absence guard."""
import pytest

from vcf2report import config
from vcf2report.annotate import gnomad, gnomad_parquet
from vcf2report.models import Variant

duckdb = pytest.importorskip("duckdb")


def _make_parquet(tmp_path):
    p = tmp_path / "g.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 0.30, 0.50, 300, 1000, 20, 0.48, 'nfe'),
        ('chr1', 200, 'C', 'G', 0.001, 0.002, 1, 1000, 0, 0.0005, 'afr'))
        AS t(chrom, pos, ref, alt, af, af_grpmax, ac, an, nhomalt, faf95, grpmax_pop))
        TO '{p}' (FORMAT PARQUET)""")
    con.close()
    return p


@pytest.fixture
def parquet(tmp_path, monkeypatch):
    p = _make_parquet(tmp_path)
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(p))
    gnomad_parquet._reset_for_tests()
    yield p
    gnomad_parquet._reset_for_tests()


def _v(chrom, pos, ref, alt):
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt)


def test_prime_and_match(parquet):
    n = gnomad_parquet.prime([_v("1", 100, "A", "T"), _v("chr1", 200, "C", "G")])
    assert n == 2
    r = gnomad_parquet.get("1-100-A-T")
    # af is the grpmax (popmax) value the ACMG engine cites, not the overall af.
    assert r["af"] == 0.50 and r["faf95"] == 0.48 and r["pop"] == "nfe" and r["hom"] == 20


def test_covered_contig_absence_is_zero(parquet):
    gnomad_parquet.prime([_v("1", 999, "A", "T")])   # chr1 covered, no such variant
    assert gnomad_parquet.get("1-999-A-T") == {
        "af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}


def test_uncovered_contig_left_unprimed(parquet):
    # chr2 is not in the parquet, so the variant is NOT primed -> caller falls back,
    # never a fabricated absence.
    gnomad_parquet.prime([_v("2", 500, "A", "T")])
    assert gnomad_parquet.get("2-500-A-T") is None


def test_chrom_prefix_tolerated(parquet):
    gnomad_parquet.prime([_v("chr1", 100, "A", "T")])
    assert gnomad_parquet.get("1-100-A-T")["af"] == 0.50


def test_off_when_unconfigured(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "GNOMAD_PARQUET", None)
    gnomad_parquet._reset_for_tests()
    assert gnomad_parquet.available() is False
    assert gnomad_parquet.prime([_v("1", 100, "A", "T")]) == 0
    gnomad_parquet._reset_for_tests()


def test_lookup_prefers_primed_parquet(parquet, monkeypatch):
    from vcf2report.annotate import cache
    monkeypatch.setattr(cache, "get", lambda *a, **k: None)
    gnomad_parquet.prime([_v("1", 100, "A", "T")])
    r = gnomad.lookup(_v("1", 100, "A", "T"))
    assert r["af"] == 0.50 and "parquet" in r["_source"]
