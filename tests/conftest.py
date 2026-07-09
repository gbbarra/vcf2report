"""Shared test fixtures.

Tests run **offline by default** (no network) and with an **isolated on-disk
cache** so they are hermetic and fast. Tests that exercise the live clients opt
back in with ``monkeypatch.setenv("OFFLINE", "")`` and mock the HTTP layer.
"""
import pytest

from vcf2report import config


@pytest.fixture(autouse=True)
def hermetic_env(tmp_path, monkeypatch):
    monkeypatch.setenv("OFFLINE", "1")
    # Default to the pure Python reader for deterministic, header-tolerant parsing;
    # the cyvcf2 path has its own dedicated tests (test_cyvcf2.py).
    monkeypatch.setenv("VCF2REPORT_NO_CYVCF2", "1")
    monkeypatch.setattr(config, "CACHE_DIR", tmp_path / "cache")
    # Never touch a developer's local ~1 GB AlphaMissense file: keep tests hermetic
    # and fast (the client degrades to "no score", as on a fresh checkout). Tests
    # that exercise AlphaMissense inject scores directly or mock the client.
    monkeypatch.setattr(config, "ALPHAMISSENSE_LOCAL", tmp_path / "no_alphamissense.tsv.gz")
    from vcf2report.annotate import alphamissense
    monkeypatch.setattr(alphamissense, "_tabix", None)
    monkeypatch.setattr(alphamissense, "_tabix_tried", False)
    yield
