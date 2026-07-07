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
    monkeypatch.setattr(config, "CACHE_DIR", tmp_path / "cache")
    yield
