"""On-disk annotation cache.

Every annotator checks the cache first. ``warm_cache.py`` pre-fills it so a live
demo is network-independent while the code remains *capable* of real calls.
Set ``OFFLINE=1`` to force cache/local-only (no network at all).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from ..config import CACHE_DIR


def _cache_file(source: str) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    return CACHE_DIR / f"{source}.json"


def load(source: str) -> dict[str, Any]:
    fp = _cache_file(source)
    if fp.exists():
        try:
            return json.loads(fp.read_text())
        except (json.JSONDecodeError, OSError):
            return {}
    return {}


def get(source: str, key: str) -> Optional[dict]:
    return load(source).get(key)


def put(source: str, key: str, value: dict) -> None:
    data = load(source)
    data[key] = value
    _cache_file(source).write_text(json.dumps(data, indent=2, sort_keys=True))
