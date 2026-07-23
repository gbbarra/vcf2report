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

# gnomAD missense-constraint thresholds.
#   PP2 fires when a gene is significantly depleted of missense variation
#   (missense z-score >= 3.09, the ClinGen SVI cut used across variant tools).
#   BP1's proxy for "primarily truncating variants cause disease" is a gene that
#   is LoF-intolerant yet shows NO missense depletion — its observed/expected
#   missense upper CI sits at/above 1.0 (missense is tolerated).
MIS_Z_CONSTRAINED = 3.09
OE_MIS_TOLERANT = 1.0


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
                mis_z = float(parts[3]) if len(parts) > 3 and parts[3] else None
                oe_mis_upper = float(parts[4]) if len(parts) > 4 and parts[4] else None
                # LoF-intolerant per gnomAD convention: pLI>=0.9 or LOEUF<0.35.
                lof_intolerant = (pli is not None and pli >= 0.9) or (
                    loeuf is not None and loeuf < 0.35)
                # Missense-constrained (PP2): significantly depleted of missense.
                missense_constrained = mis_z is not None and mis_z >= MIS_Z_CONSTRAINED
                # Missense-tolerant (BP1): no missense depletion (obs/exp CI >= 1).
                missense_tolerant = oe_mis_upper is not None and oe_mis_upper >= OE_MIS_TOLERANT
                d[gene] = {"pli": pli, "loeuf": loeuf, "lof_intolerant": lof_intolerant,
                           "mis_z": mis_z, "oe_mis_upper": oe_mis_upper,
                           "missense_constrained": missense_constrained,
                           "missense_tolerant": missense_tolerant}
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


_CONSTRAINT_NULL = {"pli": None, "loeuf": None, "lof_intolerant": None,
                    "mis_z": None, "oe_mis_upper": None,
                    "missense_constrained": None, "missense_tolerant": None}


def gene_constraint(gene: Optional[str]) -> dict:
    if not gene:
        return {**_CONSTRAINT_NULL, "_source": "gnomAD constraint (no gene)"}
    row = _load_constraint().get(gene)
    if row is None:
        return {**_CONSTRAINT_NULL, "_source": "gnomAD constraint (gene not found)"}
    return {**row, "_source": "gnomAD v2.1.1 constraint (local)"}


def insilico(variant: Variant) -> dict:
    row = _load_insilico().get(variant.key)
    if row is None:
        return {"revel": None, "cadd": None, "_source": "in-silico (none)"}
    return {**row, "_source": "REVEL/CADD (local)"}
