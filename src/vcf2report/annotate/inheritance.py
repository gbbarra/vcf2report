"""Gene -> mode of inheritance, derived from the local HPO gene-to-phenotype table.

The HPO annotates each disease gene with its inheritance mode as an ordinary phenotype
term (HP:0000006 autosomal dominant, HP:0000007 autosomal recessive, HP:0001417 X-linked
...), so ``genes_to_phenotype.tsv`` — already installed for PP4 — carries inheritance for
~4.7k genes offline, at no extra download.

This matters twice:

* **Frequency thresholds (PM2/BS1).** The credible AF ceiling is disorder-dependent: a
  recessive condition tolerates a far higher AF than a dominant one, because carriers are
  common and healthy. Without inheritance every gene falls back to the strict default.
* **PVS1.** ClinGen SVI asks whether *LoF is a known mechanism of disease for this gene*.
  The engine's original proxy for that was population constraint (pLI / LOEUF) — which
  measures selection against **heterozygous** LoF. In a recessive disorder the carrier is
  healthy, so there is no heterozygous selection and the gene scores as unconstrained
  even though LoF is precisely the mechanism. Gating PVS1 on constraint therefore blocks
  recessive disease genes by construction. (It also misfires on late-onset or
  incompletely penetrant dominants: TP53's LOEUF is 0.469, i.e. "not LoF-intolerant".)

**This is a proxy, not curation.** "The gene has an established autosomal-recessive
phenotype" is strong evidence that LoF causes disease there, but it is not the same as
ClinGen's gene-disease-validity + dosage-sensitivity curation, which is what a clinical
deployment should key PVS1 on. It replaces a proxy that is systematically wrong for
recessive genes with one that is usually right; both are labelled as such in the report.
"""
from __future__ import annotations

from typing import Optional

from .. import config
from .hpo import _read_lines

# HPO inheritance terms -> the engine's moi vocabulary. X-linked recessive/dominant both
# collapse to "XL": the engine only distinguishes recessive-style carrier frequencies from
# dominant-style near-absence, and a hemizygous male is dominant-like either way.
_TERMS = {
    "HP:0000006": "AD",    # autosomal dominant inheritance
    "HP:0000007": "AR",    # autosomal recessive inheritance
    "HP:0001417": "XL",    # X-linked inheritance
    "HP:0001419": "XL",    # X-linked recessive inheritance
    "HP:0001423": "XL",    # X-linked dominant inheritance
}

_gene_moi: Optional[dict[str, frozenset]] = None


def _load() -> dict[str, frozenset]:
    global _gene_moi
    if _gene_moi is None:
        d: dict[str, set[str]] = {}
        fp = config.HPO_GENES_LOCAL
        if fp.exists():
            for line in _read_lines(fp):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                moi = _TERMS.get(parts[1])
                if moi:
                    d.setdefault(parts[0].upper(), set()).add(moi)
        _gene_moi = {g: frozenset(v) for g, v in d.items()}  # publish only when fully built
    return _gene_moi


def modes(gene: Optional[str]) -> frozenset:
    """Every inheritance mode HPO records for the gene ({"AD","AR"} for genes with both)."""
    if not gene:
        return frozenset()
    return _load().get(gene.upper(), frozenset())


def label(gene: Optional[str]) -> Optional[str]:
    """Human-readable inheritance for the report: "AD", "AR", "AD+AR", ... or None.

    For display and audit only — never for a threshold. A gene with two inheritance modes
    must NOT be collapsed to one value before a criterion picks it up: "conservative" has a
    DIRECTION, and the two directions disagree. See config.pm2_af_ceiling / bs1_af_cutoff,
    which each take their own pole from modes().
    """
    m = modes(gene)
    return "+".join(sorted(m)) if m else None


def lof_is_disease_mechanism(gene: Optional[str]) -> bool:
    """True when HPO records an autosomal-recessive phenotype for the gene.

    Used to open the PVS1 gate for recessive genes, which population constraint cannot
    see (see the module docstring). Deliberately NOT extended to dominant genes: a
    dominant phenotype can be driven by gain-of-function or a dominant-negative just as
    easily as by haploinsufficiency, so for those the constraint evidence still has to
    carry PVS1. Recessive disease, by contrast, is overwhelmingly loss-of-function.

    ACMG classifies the VARIANT, not the patient: a heterozygous null in a recessive gene
    is still a pathogenic variant — the carrier merely needs a second hit to be affected.
    """
    return "AR" in modes(gene)
