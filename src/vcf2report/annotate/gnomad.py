"""gnomAD population-frequency client.

Resolution order: memory cache -> on-disk cache -> live GraphQL API (if online
and httpx available) -> bundled local snapshot. Returns popmax AF, AC/AN, and
homozygote count. All fields the report/ACMG engine cite come from here.
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
    genome { af ac an homozygote_count populations { id af ac an } }
    exome  { af ac an homozygote_count populations { id af ac an } }
  }
}
"""


def _load_local() -> dict:
    global _local
    if _local is None:
        fp = config.GNOMAD_LOCAL
        _local = json.loads(fp.read_text()) if fp.exists() else {}
    return _local


def _from_payload(payload: dict) -> dict:
    """Reduce a gnomAD variant payload to popmax AF + counts."""
    best = {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None}
    for src in ("exome", "genome"):
        block = payload.get(src)
        if not block:
            continue
        best["hom"] = max(best["hom"], block.get("homozygote_count") or 0)
        for pop in block.get("populations", []) or []:
            if (pop.get("af") or 0.0) > best["af"]:
                best = {"af": pop["af"], "ac": pop.get("ac", 0),
                        "an": pop.get("an", 0), "hom": best["hom"], "pop": pop.get("id")}
    return best


def _live(variant: Variant) -> Optional[dict]:  # pragma: no cover - network
    try:
        import httpx
    except ImportError:
        return None
    try:
        resp = httpx.post(
            config.GNOMAD_API,
            json={"query": _GRAPHQL,
                  "variables": {"variantId": variant.key, "dataset": config.GNOMAD_DATASET}},
            timeout=15.0,
        )
        resp.raise_for_status()
        data = resp.json().get("data", {}).get("variant")
        if not data:
            return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None}
        return _from_payload(data)
    except Exception:
        return None


def lookup(variant: Variant) -> dict:
    """Return {'af','ac','an','hom','pop','_source'} for a variant."""
    cached = cache.get(_SOURCE, variant.key)
    if cached is not None:
        return {**cached, "_source": f"gnomAD {config.GNOMAD_DATASET} (cache)"}

    if not config.offline():
        live = _live(variant)
        if live is not None:
            cache.put(_SOURCE, variant.key, live)
            return {**live, "_source": f"gnomAD {config.GNOMAD_DATASET} (live)"}

    local = _load_local().get(variant.key)
    if local is not None:
        return {**local, "_source": f"gnomAD {config.GNOMAD_DATASET} (local snapshot)"}

    # Not found anywhere == absent from gnomAD.
    return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None,
            "_source": f"gnomAD {config.GNOMAD_DATASET} (not observed)"}
