"""ABraOM (Brazilian) allele-frequency lookup.

ABraOM is the SABE cohort of admixed Brazilian genomes (http://abraom.ib.usp.br).
There is no live API — it ships as a static dataset. Consulting it alongside
gnomAD is the local differentiator: a variant absent from gnomAD but common in
Brazilians is correctly down-weighted (blocks PM2, can trigger BA1/BS1),
preventing a real class of misclassifications in admixed patients.
"""
from __future__ import annotations

from typing import Optional

from .. import config
from ..models import Variant

_local: Optional[dict] = None
_COLUMNS = ["key", "af", "ac", "an"]


def _load_local() -> dict:
    global _local
    if _local is None:
        _local = {}
        fp = config.ABRAOM_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                row = dict(zip(_COLUMNS, parts))
                if row.get("key"):
                    _local[row["key"]] = row
    return _local


def lookup(variant: Variant) -> dict:
    row = _load_local().get(variant.key)
    if row is not None:
        return {"af": float(row.get("af", 0.0)), "ac": int(row.get("ac", 0) or 0),
                "an": int(row.get("an", 0) or 0), "_source": "ABraOM SABE (local)"}
    return {"af": 0.0, "ac": 0, "an": 0, "_source": "ABraOM SABE (not observed)"}
