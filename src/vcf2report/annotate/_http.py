"""Tiny HTTP helper for the live annotation clients.

Uses only the standard library (``urllib``) so live API calls work with zero
extra dependencies. Handles JSON GET/POST with timeouts, polite retry/backoff on
429/5xx, and a shared minimum-interval throttle (NCBI wants <=3 req/s, or 10 with
an API key). All errors are swallowed into ``None`` so callers fall back cleanly.
"""
from __future__ import annotations

import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Optional

_USER_AGENT = "vcf2report/0.1 (+https://github.com/gbbarra/vcf2report)"
_last_call: dict[str, float] = {}


def throttle(key: str, min_interval: float) -> None:
    """Block until at least ``min_interval`` seconds have passed for ``key``."""
    now = time.monotonic()
    prev = _last_call.get(key)
    if prev is not None:
        wait = min_interval - (now - prev)
        if wait > 0:
            time.sleep(wait)
    _last_call[key] = time.monotonic()


def _request(req: urllib.request.Request, timeout: float, retries: int) -> Optional[dict]:
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read().decode("utf-8")
            return json.loads(raw)
        except urllib.error.HTTPError as e:
            # Retry on rate-limit / transient server errors with backoff.
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return None
        except (urllib.error.URLError, TimeoutError, ValueError, OSError):
            if attempt < retries:
                time.sleep(2 ** attempt)
                continue
            return None
    return None


def post_json(url: str, payload: dict, timeout: float = 15.0,
              retries: int = 2) -> Optional[dict]:
    data = json.dumps(payload).encode("utf-8")
    req = urllib.request.Request(
        url, data=data, method="POST",
        headers={"Content-Type": "application/json", "Accept": "application/json",
                 "User-Agent": _USER_AGENT},
    )
    return _request(req, timeout, retries)


def get_json(url: str, params: dict, timeout: float = 15.0,
             retries: int = 2) -> Optional[dict]:
    qs = urllib.parse.urlencode({k: v for k, v in params.items() if v is not None})
    req = urllib.request.Request(
        f"{url}?{qs}", method="GET",
        headers={"Accept": "application/json", "User-Agent": _USER_AGENT},
    )
    return _request(req, timeout, retries)
