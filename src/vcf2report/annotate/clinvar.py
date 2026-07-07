"""ClinVar clinical-significance client.

Primary source is a bundled ClinVar slice (TSV keyed by CHROM-POS-REF-ALT) —
offline, deterministic, fast. Live NCBI E-utilities is used as a flourish when
online and httpx is available. Returns significance, review status, accession,
condition, and date.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

from .. import config
from ..models import Variant
from . import cache

_SOURCE = "clinvar"
_local: Optional[dict] = None
_COLUMNS = ["key", "significance", "review_status", "accession", "condition", "date"]


def _load_local() -> dict:
    global _local
    if _local is None:
        _local = {}
        fp = config.CLINVAR_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                row = dict(zip(_COLUMNS, parts))
                if row.get("key"):
                    _local[row["key"]] = row
    return _local


def _live(variant: Variant) -> Optional[dict]:  # pragma: no cover - network
    try:
        import httpx  # noqa: F401
    except ImportError:
        return None
    # E-utilities esearch+esummary would go here (db=clinvar). Intentionally a
    # thin stub: the bundled slice is authoritative for the demo and E-utilities
    # is rate-limited (3/s; 10/s with NCBI_API_KEY). Returning None falls back.
    return None


def lookup(variant: Variant) -> dict:
    cached = cache.get(_SOURCE, variant.key)
    if cached is not None:
        return {**cached, "_source": "ClinVar (cache)"}

    if not config.offline():
        live = _live(variant)
        if live is not None:
            cache.put(_SOURCE, variant.key, live)
            return {**live, "_source": "ClinVar (live E-utilities)"}

    local = _load_local().get(variant.key)
    if local is not None:
        date = local.get("date", "")
        return {**local, "_source": f"ClinVar slice ({date})"}

    return {"significance": None, "review_status": None, "accession": None,
            "condition": None, "date": None, "_source": "ClinVar (no record)"}
