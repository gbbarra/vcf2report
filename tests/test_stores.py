"""Store health / integrity / freshness — the checker must actually catch problems:
a missing store, a corrupt (unreadable) store, an incomplete one (missing a core contig or a
row count that disagrees with the build manifest), and a stale weekly store past its window.
"""
from datetime import datetime, timedelta, timezone

import pytest

duckdb = pytest.importorskip("duckdb")

from vcf2report import config, stores


def _mk_clinvar(tmp_path, chroms, rows_per=1):
    """A ClinVar-shaped partitioned parquet with the given contigs."""
    p = tmp_path / "cv_parquet"
    rows = []
    for c in chroms:
        for i in range(rows_per):
            rows.append(f"('{c}', {100 + i}, 'A', 'T', 'Pathogenic', "
                        f"'reviewed by expert panel', 3, 'VCV{i}', 'Cond')")
    con = duckdb.connect()
    con.execute(f"""COPY (SELECT * FROM (VALUES {','.join(rows)})
        t(chrom,pos,ref,alt,significance,review_status,review_stars,accession,condition))
        TO '{p}' (FORMAT PARQUET, PARTITION_BY (chrom))""")
    con.close()
    return p


def _use(monkeypatch, path):
    monkeypatch.setattr(config, "CLINVAR_PARQUET", path)
    stores._reset_for_tests()


def test_healthy_store_ok(tmp_path, monkeypatch):
    _use(monkeypatch, _mk_clinvar(tmp_path, sorted(stores.CORE_CHROMS)))
    stores.write_manifest("clinvar")
    e = stores.store_health("clinvar")["clinvar"]
    assert e["status"] == "ok" and e["present"] and e["readable"] and e["complete"]
    assert e["rows"] == len(stores.CORE_CHROMS) and not e["update_recommended"]


def test_missing_store(tmp_path, monkeypatch):
    _use(monkeypatch, tmp_path / "nope")
    e = stores.store_health("clinvar")["clinvar"]
    assert e["status"] == "missing" and not e["present"] and e["update_recommended"]


def test_incomplete_missing_core_contig(tmp_path, monkeypatch):
    chroms = sorted(stores.CORE_CHROMS - {"chrX"})   # drop a required contig
    _use(monkeypatch, _mk_clinvar(tmp_path, chroms))
    stores.write_manifest("clinvar")
    e = stores.store_health("clinvar")["clinvar"]
    assert e["status"] == "incomplete" and e["complete"] is False
    assert "chrX" in e["missing_core_chroms"] and e["update_recommended"]


def test_row_count_mismatch_flags_incomplete(tmp_path, monkeypatch):
    _use(monkeypatch, _mk_clinvar(tmp_path, sorted(stores.CORE_CHROMS)))
    stores.write_manifest("clinvar")
    # Corrupt the manifest's row count -> the measured store must disagree.
    import json
    mf = config.CLINVAR_PARQUET / stores.MANIFEST
    d = json.loads(mf.read_text())
    d["rows"] = 999999
    mf.write_text(json.dumps(d))
    e = stores.store_health("clinvar")["clinvar"]
    assert e["status"] == "incomplete" and e["complete"] is False
    assert e["rows_mismatch"]["expected"] == 999999


def test_stale_weekly_store_flagged(tmp_path, monkeypatch):
    _use(monkeypatch, _mk_clinvar(tmp_path, sorted(stores.CORE_CHROMS)))
    old = (datetime.now(timezone.utc) - timedelta(days=30)).isoformat(timespec="seconds")
    stores.write_manifest("clinvar", built_at=old)
    e = stores.store_health("clinvar")["clinvar"]
    assert e["status"] == "stale" and e["update_recommended"]
    assert e["age_days"] >= 29 and "weekly" in e["reason"]


def test_quick_mode_skips_row_scan(tmp_path, monkeypatch):
    _use(monkeypatch, _mk_clinvar(tmp_path, sorted(stores.CORE_CHROMS)))
    stores.write_manifest("clinvar")
    e = stores.store_health("clinvar", measure=False)["clinvar"]
    assert e["present"] and e.get("rows") is None and e["status"] == "ok"


def test_gate_blocks_on_missing_or_corrupt(monkeypatch):
    monkeypatch.setattr(stores, "store_health", lambda measure=True: {
        "gnomad": {"status": "ok"}, "alphamissense": {"status": "missing"},
        "clinvar": {"status": "stale"}})
    g = stores.gate()
    assert not g["ready"] and g["blocking"] == ["alphamissense"] and g["stale"] == ["clinvar"]
    monkeypatch.setattr(stores, "store_health", lambda measure=True: {
        "gnomad": {"status": "corrupt"}, "alphamissense": {"status": "ok"},
        "clinvar": {"status": "ok"}})
    assert not stores.gate()["ready"]


def test_gate_ready_when_intact_even_if_stale(monkeypatch):
    # A stale-but-intact ClinVar warns but must NOT block the analysis.
    monkeypatch.setattr(stores, "store_health", lambda measure=True: {
        "gnomad": {"status": "ok"}, "alphamissense": {"status": "ok"},
        "clinvar": {"status": "stale"}})
    g = stores.gate()
    assert g["ready"] and g["blocking"] == [] and g["stale"] == ["clinvar"]


def test_frozen_store_not_stale_without_manifest(tmp_path, monkeypatch):
    # gnomAD/AM are frozen: no manifest date -> still not flagged for a routine update.
    p = _mk_clinvar(tmp_path, sorted(stores.CORE_CHROMS))  # shape irrelevant; reuse
    monkeypatch.setattr(config, "ALPHAMISSENSE_PARQUET", p)
    stores._reset_for_tests()
    e = stores.store_health("alphamissense", measure=False)["alphamissense"]
    assert e["cadence"] == "frozen" and not e["update_recommended"]
