"""Gene-level constraint and in-silico predictor lookups (local datasets).

* Gene constraint (gnomAD-derived pLI/LOEUF) drives PVS1's "LoF-intolerant gene"
  requirement and can inform PP2/BP1.
* In-silico scores (REVEL, CADD) drive PP3/BP4.

Both ship as small local TSVs for a deterministic, offline demo.
"""
from __future__ import annotations

from typing import Optional

from .. import config
from ..models import Variant

_constraint: Optional[dict] = None
_insilico: Optional[dict] = None


def _read_lines(fp):
    """Yield lines from a plain or gzip-compressed text file."""
    if str(fp).endswith(".gz"):
        import gzip
        with gzip.open(fp, "rt") as fh:
            for line in fh:
                yield line.rstrip("\n")
    else:
        yield from fp.read_text().splitlines()


def _load_constraint() -> dict:
    global _constraint
    if _constraint is None:
        d: dict = {}
        fp = config.CONSTRAINT_LOCAL
        if fp.exists():
            for line in _read_lines(fp):
                if not line.strip() or line.startswith("#") or line.startswith("gene\t"):
                    continue
                parts = line.split("\t")
                gene = parts[0]
                pli = float(parts[1]) if len(parts) > 1 and parts[1] else None
                loeuf = float(parts[2]) if len(parts) > 2 and parts[2] else None
                # LoF-intolerant per gnomAD convention: pLI>=0.9 or LOEUF<0.35.
                lof_intolerant = (pli is not None and pli >= 0.9) or (
                    loeuf is not None and loeuf < 0.35)
                d[gene] = {"pli": pli, "loeuf": loeuf, "lof_intolerant": lof_intolerant}
        _constraint = d  # publish only when fully built
    return _constraint


def _load_insilico() -> dict:
    global _insilico
    if _insilico is None:
        d: dict = {}
        fp = config.INSILICO_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#") or line.startswith("key\t"):
                    continue
                parts = line.split("\t")
                key = parts[0]
                d[key] = {
                    "revel": float(parts[1]) if len(parts) > 1 and parts[1] else None,
                    "cadd": float(parts[2]) if len(parts) > 2 and parts[2] else None,
                }
        _insilico = d  # publish only when fully built
    return _insilico


def gene_constraint(gene: Optional[str]) -> dict:
    if not gene:
        return {"pli": None, "loeuf": None, "lof_intolerant": None,
                "_source": "gnomAD constraint (no gene)"}
    row = _load_constraint().get(gene)
    if row is None:
        return {"pli": None, "loeuf": None, "lof_intolerant": None,
                "_source": "gnomAD constraint (gene not found)"}
    return {**row, "_source": "gnomAD v2.1.1 LoF constraint (local)"}


def insilico(variant: Variant) -> dict:
    row = _load_insilico().get(variant.key)
    if row is None:
        return {"revel": None, "cadd": None, "_source": "in-silico (none)"}
    return {**row, "_source": "REVEL/CADD (local)"}
