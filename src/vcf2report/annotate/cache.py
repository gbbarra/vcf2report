"""On-disk annotation cache.

Every annotator checks the cache first. ``warm_cache.py`` pre-fills it so a live
demo is network-independent while the code remains *capable* of real calls.
Set ``OFFLINE=1`` to force cache/local-only (no network at all).
"""
from __future__ import annotations

import json
import os
import tempfile
import threading
from pathlib import Path
from typing import Any, Optional

from .. import config

_write_lock = threading.Lock()


def _cache_file(source: str) -> Path:
    # Read config.CACHE_DIR at call time so overrides (env / tests) take effect.
    config.CACHE_DIR.mkdir(parents=True, exist_ok=True)
    # The cache holds patient-derived variant coordinates — keep it owner-only.
    try:
        os.chmod(config.CACHE_DIR, 0o700)
    except OSError:
        pass
    return config.CACHE_DIR / f"{source}.json"


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
    # Lock the whole read-modify-write so concurrent writers can't lose entries
    # (last-writer-wins). Atomic replace also prevents a truncated file.
    with _write_lock:
        data = load(source)
        data[key] = value
        fp = _cache_file(source)
        fd, tmp = tempfile.mkstemp(dir=str(fp.parent), suffix=".tmp")
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(data, fh, indent=2, sort_keys=True)
            os.replace(tmp, fp)
        except OSError:
            if os.path.exists(tmp):
                os.unlink(tmp)
