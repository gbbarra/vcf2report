"""gnomAD population-frequency client.

Resolution order: on-disk cache -> live GraphQL API (unless OFFLINE) -> bundled
local snapshot. Returns popmax AF, AC/AN, and homozygote count. All fields the
report/ACMG engine cite come from here.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from .. import config
from ..models import Variant
from . import cache

_SOURCE = "gnomad"
_local: Optional[dict] = None

_GRAPHQL = """
query Variant($variantId: String!, $dataset: DatasetId!) {
  variant(variantId: $variantId, dataset: $dataset) {
    variant_id
    genome { ac an homozygote_count populations { id ac an } }
    exome  { ac an homozygote_count populations { id ac an } }
  }
}
"""


def _load_local() -> dict:
    global _local
    if _local is None:
        fp = config.GNOMAD_LOCAL
        _local = json.loads(fp.read_text()) if fp.exists() else {}
    return _local


# Population ids to exclude from popmax (bottleneck / non-continental groups),
# Groups gnomAD excludes from grpmax/popmax (bottleneck / small cohorts). The
# large continental ancestry groups — afr, amr, eas, nfe, sas — are all INCLUDED.
_POPMAX_EXCLUDE = {"asj", "fin", "oth", "remaining", "mid", "ami"}
# Ignore degenerate tiny cohorts when picking popmax (guards 1/2-allele artefacts;
# well below any real continental group's sample size).
_MIN_AN = 100


def _af(ac: int, an: int) -> float:
    return (ac / an) if an else 0.0


def _from_payload(payload: dict) -> dict:
    """Reduce a gnomAD variant payload to popmax AF + counts.

    AC/AN are summed per population ACROSS exome+genome (the callsets are largely
    non-overlapping) before computing AF, and homozygote counts are summed. popmax
    is the highest per-population AF over the included continental groups.
    """
    pop_ac: dict[str, int] = {}
    pop_an: dict[str, int] = {}
    total_hom = 0
    for src in ("exome", "genome"):
        block = payload.get(src)
        if not block:
            continue
        total_hom += block.get("homozygote_count") or 0
        for pop in block.get("populations", []) or []:
            pid = (pop.get("id") or "").lower()
            # Skip sex-stratified / sub-population ids (e.g. "nfe_XX") and the
            # excluded bottleneck groups; keep top-level ancestry groups only.
            if "_" in pid or pid in _POPMAX_EXCLUDE:
                continue
            pop_ac[pid] = pop_ac.get(pid, 0) + (pop.get("ac") or 0)
            pop_an[pid] = pop_an.get(pid, 0) + (pop.get("an") or 0)

    best = {"af": 0.0, "ac": 0, "an": 0, "hom": total_hom, "pop": None}
    for pid, an in pop_an.items():
        if an < _MIN_AN:
            continue
        af = _af(pop_ac[pid], an)
        if af > best["af"]:
            best = {"af": af, "ac": pop_ac[pid], "an": an, "hom": total_hom, "pop": pid}
    return best


def _is_not_found(errors: list) -> bool:
    """True only if every GraphQL error is the benign 'Variant not found'."""
    return bool(errors) and all(
        "not found" in (e.get("message") or "").lower() for e in errors
    )


def _live(variant: Variant) -> Optional[dict]:
    """Query the gnomAD GraphQL API.

    Returns an absent-record dict only for a genuine 'variant not found' (clean
    response). A real failure — transport error, or a 200 carrying a non-"not
    found" GraphQL ``errors`` block — returns None so the caller falls back to the
    local snapshot instead of caching a fabricated AF 0.
    """
    from . import _http

    resp = _http.post_json(
        config.GNOMAD_API,
        {"query": _GRAPHQL,
         "variables": {"variantId": variant.key, "dataset": config.GNOMAD_DATASET}},
        timeout=15.0,
    )
    if resp is None:
        return None  # transport error -> fall back
    errors = resp.get("errors") or []
    if errors and not _is_not_found(errors):
        return None  # server-side failure (rate limit, timeout, schema) -> fall back
    if "data" not in resp:
        return None
    data = (resp.get("data") or {}).get("variant")
    if not data:
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None}  # genuinely absent
    return _from_payload(data)


def _unknown(reason: str) -> dict:
    """Sentinel for 'frequency unavailable' — distinct from a confirmed absence.

    AF is None so downstream ACMG logic does NOT treat it as rare/absent (PM2 must
    not fire on unknown frequency).
    """
    return {"af": None, "ac": None, "an": None, "hom": None, "pop": None,
            "_source": f"gnomAD {config.GNOMAD_DATASET} ({reason})"}


def lookup(variant: Variant) -> dict:
    """Return {'af','ac','an','hom','pop','_source'} for a variant."""
    cached = cache.get(_SOURCE, variant.key)
    if cached is not None:
        return {**cached, "_source": f"gnomAD {config.GNOMAD_DATASET} (cache)"}

    live_failed = False
    if not config.offline():
        live = _live(variant)
        if live is not None:
            cache.put(_SOURCE, variant.key, live)
            return {**live, "_source": f"gnomAD {config.GNOMAD_DATASET} (live)"}
        live_failed = True  # transport/server failure -> try local, else unknown

    local = _load_local().get(variant.key)
    if local is not None:
        return {**local, "_source": f"gnomAD {config.GNOMAD_DATASET} (local snapshot)"}

    # No definitive answer: a failed live lookup or offline with no local record.
    # Report 'unavailable' (AF None), NOT a fabricated absence.
    return _unknown("unavailable — lookup failed" if live_failed else "no local data")
