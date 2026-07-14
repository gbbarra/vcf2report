"""AlphaMissense + ClinVar Parquet stores wired into the annotate stage.

Builds tiny single-file parquets in a tmp dir and checks: the batch prime resolves hits,
never asserts AM absence (miss -> None), a ClinVar miss falls through (get -> None), the
clients prefer the parquet, and everything degrades to the tabix path when the store /
duckdb is absent.
"""
import pytest

duckdb = pytest.importorskip("duckdb")

from vcf2report import config
from vcf2report.annotate import (alphamissense, alphamissense_parquet, clinvar,
                                 clinvar_parquet)
from vcf2report.annotate import cache as anncache
from vcf2report.models import Variant


def _write_am(tmp_path):
    p = tmp_path / "am.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 0.9::DOUBLE, 'likely_pathogenic'),
        ('chr1', 200, 'C', 'G', 0.1::DOUBLE, 'likely_benign')
      ) t(chrom,pos,ref,alt,am_pathogenicity,am_class)) TO '{p}' (FORMAT PARQUET)""")
    con.close()
    return p


def _write_cv(tmp_path):
    p = tmp_path / "cv.parquet"
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES
        ('chr1', 100, 'A', 'T', 'Pathogenic', 'reviewed by expert panel', 3, 'VCV1', 'Cond')
      ) t(chrom,pos,ref,alt,significance,review_status,review_stars,accession,condition))
      TO '{p}' (FORMAT PARQUET)""")
    con.close()
    return p


def test_am_parquet_prime_hit_and_absence(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ALPHAMISSENSE_PARQUET", _write_am(tmp_path))
    alphamissense_parquet._reset_for_tests()
    assert alphamissense_parquet.available()
    hit = Variant(chrom="1", pos=100, ref="A", alt="T")
    miss = Variant(chrom="1", pos=999, ref="A", alt="T")
    got = alphamissense_parquet.prime([hit, miss])
    assert got[hit.key] == {"am_pathogenicity": 0.9, "am_class": "likely_pathogenic"}
    assert got[miss.key] is None  # AM absence is never asserted -> None (fires nothing)


def test_alphamissense_prime_prefers_parquet(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ALPHAMISSENSE_PARQUET", _write_am(tmp_path))
    alphamissense_parquet._reset_for_tests()
    alphamissense._reset_for_tests()
    v = Variant(chrom="1", pos=100, ref="A", alt="T")
    assert alphamissense.prime([v]) == 1
    r = alphamissense.lookup(v)
    assert r["am_pathogenicity"] == 0.9 and "primed" in r["_source"]


def test_clinvar_parquet_prime_and_get(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "CLINVAR_PARQUET", _write_cv(tmp_path))
    clinvar_parquet._reset_for_tests()
    v = Variant(chrom="1", pos=100, ref="A", alt="T")
    assert clinvar_parquet.prime([v]) == 1
    rec = clinvar_parquet.get(v.key)
    assert rec["significance"] == "Pathogenic" and rec["review_status"] == "reviewed by expert panel"
    assert clinvar_parquet.get("9-9-A-T") is None  # miss -> None -> lookup falls through


def test_clinvar_lookup_uses_parquet(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "CLINVAR_PARQUET", _write_cv(tmp_path))
    monkeypatch.setattr(anncache, "get", lambda *a, **k: None)
    clinvar_parquet._reset_for_tests()
    v = Variant(chrom="1", pos=100, ref="A", alt="T")
    clinvar_parquet.prime([v])
    r = clinvar.lookup(v)
    assert r["significance"] == "Pathogenic" and r["_source"] == "ClinVar (local)"


def test_parquet_unavailable_falls_back(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "ALPHAMISSENSE_PARQUET", tmp_path / "nope")
    monkeypatch.setattr(config, "CLINVAR_PARQUET", tmp_path / "nope2")
    alphamissense_parquet._reset_for_tests()
    clinvar_parquet._reset_for_tests()
    v = Variant(chrom="1", pos=1, ref="A", alt="T")
    assert not alphamissense_parquet.available() and alphamissense_parquet.prime([v]) is None
    assert not clinvar_parquet.available() and clinvar_parquet.prime([v]) == 0
