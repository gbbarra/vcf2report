"""Allele frequency from a LOCAL cohort the operator supplies.

Population databases are dominated by a few ancestries, so a variant that is rare or
absent in them can still be common in the population a lab actually serves — and a rare
allele is exactly what earns PM2. Consulting a local cohort alongside gnomAD catches that
class of misclassification: a locally common variant is blocked from PM2 and can trigger
BA1/BS1.

**This project ships no cohort.** The mechanism is here; the data is yours to supply,
under whatever terms govern it. Nothing is downloaded and nothing is redistributed.

Two ways to provide it:

* a TSV at ``VCF2REPORT_LOCAL_COHORT`` (default ``data/local_cohort/frequencies.tsv``),
  tab-separated ``key<TAB>af<TAB>ac<TAB>an`` where ``key`` is ``chrom-pos-ref-alt`` with
  no ``chr`` prefix and upper-case alleles (matching :attr:`Variant.key`);
* a ``LOCAL_AF`` INFO field written by your own annotation step (see
  ``config.INFO_ALIASES``).

A variant absent from the table is reported as **unknown**, never as a checked zero —
the table is a slice of one cohort, so a miss means "not in this table", not "confirmed
absent locally". The criteria that read it say which databases were actually consulted.
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
        d: dict = {}
        fp = config.LOCAL_COHORT
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                row = dict(zip(_COLUMNS, parts))
                if row.get("key"):
                    d[row["key"]] = row
        _local = d  # publish only when fully built
    return _local


def _num(row: dict, field: str, cast):
    """A missing or blank cell is UNKNOWN, not zero.

    A truncated row previously became a surveyed zero stamped '(local)' — the frequency
    miss this feature exists to prevent — and an empty AF cell raised ValueError out of
    the annotation pass.
    """
    raw = (row.get(field) or "").strip()
    if not raw:
        return None
    try:
        return cast(raw)
    except ValueError:
        return None


def lookup(variant: Variant) -> dict:
    row = _load_local().get(variant.key)
    if row is not None:
        af = _num(row, "af", float)
        if af is None:                      # a row with no usable AF answers nothing
            return {"af": None, "ac": None, "an": None,
                    "_source": "local cohort (row present, no usable AF)"}
        return {"af": af, "ac": _num(row, "ac", int), "an": _num(row, "an", int),
                "_source": "local cohort (local table)"}
    return {"af": None, "ac": None, "an": None,
            "_source": "local cohort (not in local table)"}
