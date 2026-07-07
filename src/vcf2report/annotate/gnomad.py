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
# mirroring gnomAD's own popmax convention.
_POPMAX_EXCLUDE = {"asj", "fin", "oth", "remaining", "mid", "ami", "sas"}


def _af(ac: int, an: int) -> float:
    return (ac / an) if an else 0.0


def _from_payload(payload: dict) -> dict:
    """Reduce a gnomAD variant payload to popmax AF + counts.

    AF is computed from ac/an (robust to schema differences) and popmax is the
    highest per-population AF across exome+genome, excluding bottleneck groups.
    """
    best = {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None}
    for src in ("exome", "genome"):
        block = payload.get(src)
        if not block:
            continue
        best["hom"] = max(best["hom"], block.get("homozygote_count") or 0)
        for pop in block.get("populations", []) or []:
            pid = (pop.get("id") or "").lower()
            # gnomAD population ids can be suffixed (e.g. "nfe_bgr"); keep the
            # top-level ancestry groups only (no underscore) for popmax.
            if "_" in pid or pid in _POPMAX_EXCLUDE:
                continue
            ac, an = pop.get("ac") or 0, pop.get("an") or 0
            af = _af(ac, an)
            if af > best["af"]:
                best = {"af": af, "ac": ac, "an": an, "hom": best["hom"], "pop": pid}
    return best


def _live(variant: Variant) -> Optional[dict]:
    """Query the gnomAD GraphQL API. Returns absent-record dict if not found."""
    from . import _http

    resp = _http.post_json(
        config.GNOMAD_API,
        {"query": _GRAPHQL,
         "variables": {"variantId": variant.key, "dataset": config.GNOMAD_DATASET}},
        timeout=15.0,
    )
    if resp is None:
        return None  # network/transport error -> let caller fall back
    # A "variant not found" comes back as data.variant == null (often with an
    # errors block); that means the allele is absent from gnomAD, not a failure.
    data = (resp.get("data") or {}).get("variant")
    if not data:
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "pop": None}
    return _from_payload(data)


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
