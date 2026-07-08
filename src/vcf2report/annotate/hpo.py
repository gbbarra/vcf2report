"""HPO gene<->phenotype matching.

Primary source is the bundled ``genes_to_phenotype.tsv`` (from the HPO project;
columns: gene, hpo_id, hpo_name). Given a gene and the patient's HPO terms, we
compute an overlap score and the matched terms — this drives PP4 and the
candidate phenotype ranking. Live ontology.jax.org lookups can augment this but
the files are the demo-safe default.
"""
from __future__ import annotations

from typing import Optional

from .. import config

_gene_terms: Optional[dict[str, set[str]]] = None
_term_names: dict[str, str] = {}


def _load() -> dict[str, set[str]]:
    global _gene_terms
    if _gene_terms is None:
        d: dict[str, set[str]] = {}
        names: dict[str, str] = {}
        fp = config.HPO_GENES_LOCAL
        if fp.exists():
            for line in fp.read_text().splitlines():
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 2:
                    continue
                gene, hpo_id = parts[0], parts[1]
                d.setdefault(gene, set()).add(hpo_id)
                if len(parts) >= 3:
                    names[hpo_id] = parts[2]
        _term_names.update(names)
        _gene_terms = d  # publish only when fully built
    return _gene_terms


def match(gene: Optional[str], patient_hpo: list[str]) -> dict:
    """Return {'score','matched_terms','_source'} for a gene vs patient terms.

    Score is the fraction of the patient's HPO terms explained by the gene
    (|intersection| / |patient terms|), a simple, interpretable specificity proxy.
    """
    if not gene or not patient_hpo:
        return {"score": 0.0, "matched_terms": [], "_source": "HPO (no gene/terms)"}
    gene_terms = _load().get(gene, set())
    patient = set(patient_hpo)
    matched = sorted(patient & gene_terms)
    score = round(len(matched) / len(patient), 3) if patient else 0.0
    labelled = [f"{t} ({_term_names.get(t, '?')})" for t in matched]
    return {"score": score, "matched_terms": labelled,
            "_source": "HPO genes_to_phenotype (local)"}
