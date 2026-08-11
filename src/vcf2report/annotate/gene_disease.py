"""Gene -> associated diseases, with inheritance tied to each **disease** (not the gene).

For the laudo: alongside a variant's gene and ACMG tier, name the disease(s) that gene
causes and how each is inherited. Inheritance MUST be per-disease, never a single label
per gene — one gene routinely causes different diseases with different inheritance (JUP:
Naxos disease is autosomal **recessive**, but arrhythmogenic RV dysplasia 12, same gene,
is autosomal **dominant**; a merged "AD/AR" would hide exactly what matters).

Data comes from the HPO project's ``genes_to_phenotype.txt`` (gene x disease x HPO term,
where inheritance is itself an HPO term under HP:0000005) joined to ``phenotype.hpoa`` for
the readable disease name. ``scripts/build_gene_disease_store.py`` freezes that join into
``config.HPO_GENE_DISEASE_LOCAL`` (a gzipped JSON keyed by gene). No network at runtime.
"""
from __future__ import annotations

import gzip
import json
from typing import Optional

from .. import config

# HPO "Mode of inheritance" (HP:0000005) terms -> readable label for the laudo. The
# clinically-actionable mendelian modes; other subtree terms are ignored (kept out of the
# card rather than shown as noise). Shared with the store builder.
INHERITANCE_TERMS = {
    "HP:0000006": "Autosomal dominant",
    "HP:0000007": "Autosomal recessive",
    "HP:0001417": "X-linked",
    "HP:0001419": "X-linked recessive",
    "HP:0001423": "X-linked dominant",
    "HP:0001450": "Y-linked",
    "HP:0001427": "Mitochondrial",
    "HP:0010984": "Digenic",
    "HP:0032113": "Semidominant",
}

_store: Optional[dict[str, list]] = None


def _load() -> dict[str, list]:
    global _store
    if _store is None:
        fp = config.HPO_GENE_DISEASE_LOCAL
        if fp.exists():
            with gzip.open(fp, "rt") as fh:
                _store = json.load(fh)
        else:
            _store = {}
    return _store


def lookup(gene: Optional[str]) -> list[dict]:
    """Diseases for ``gene`` as ``[{disease_id, disease_name, inheritance:[...]}]``.

    Empty list when the gene is unknown or the store is not installed — the laudo simply
    omits the disease line rather than asserting the gene causes nothing.
    """
    if not gene:
        return []
    return _load().get(gene.upper(), [])
