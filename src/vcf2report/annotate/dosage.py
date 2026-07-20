"""ClinGen dosage-sensitivity — genes where loss of function is an established disease mechanism.

ClinGen SVI's PVS1 asks whether LoF is a *known mechanism of disease* for the gene. Population
constraint (pLI/LOEUF) is only a proxy for it, and a poor one for two classes:

* **recessive genes** — the carrier is healthy, so no heterozygous selection, so the gene scores
  as unconstrained even though LoF is the mechanism (handled by ``inheritance.lof_is_disease_mechanism``);
* **late-onset / incompletely penetrant dominants** — e.g. TP53 (LOEUF 0.469, "not intolerant") is a
  textbook haploinsufficient tumour suppressor, but selection is weak so constraint misses it.

ClinGen's **Haploinsufficiency = 3** ("sufficient evidence for dosage pathogenicity") is exactly the
curated gene→mechanism statement PVS1 wants — an expert panel's judgement that losing one copy causes
disease. This module exposes that set (418 genes) as a third, authoritative PVS1 mechanism route,
alongside constraint and the recessive-phenotype route.
"""
from __future__ import annotations

from typing import Optional

from .. import config

_hi: Optional[set] = None


def _load() -> set:
    global _hi
    if _hi is None:
        s: set = set()
        fp = config.CLINGEN_HI_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    s.add(line.upper())
        _hi = s  # publish only when fully built (empty set if the file is absent)
    return _hi


def haploinsufficient(gene: Optional[str]) -> bool:
    """True when ClinGen curates the gene as Haploinsufficiency=3 — LoF is an established mechanism."""
    return bool(gene) and gene.upper() in _load()
