"""HPO gene<->phenotype matching.

Primary source is the bundled ``genes_to_phenotype.tsv`` (from the HPO project;
columns: gene, hpo_id, hpo_name). Given a gene and the patient's HPO terms, we
score how well the gene's known phenotype spectrum explains the patient's terms.
This drives PP4 and the candidate phenotype ranking.

Two scorers, in order of preference:

* **ontology-aware (Lin/IC)** — when ``hpo_graph.tsv.gz`` is present (built by
  ``scripts/build_hpo_graph.py``), similarity uses the HPO ``is_a`` graph and each
  term's Information Content, so a patient term matches a *related* gene term (a
  parent/child), weighted by specificity. A best-match average over the patient's
  terms means adding more (explained) phenotypes no longer dilutes the score — the
  failure mode of exact overlap on phenotype-rich cases.
* **exact overlap** — the dependency-free fallback when the graph is absent:
  ``|patient ∩ gene_terms| / |patient|``. Unchanged behaviour, demo-safe.
"""
from __future__ import annotations

from typing import Optional

from .. import config

_gene_terms: Optional[dict[str, set[str]]] = None
_term_names: dict[str, str] = {}
# graph = (parents: id->list[id], ic: id->float, names: id->name); None until loaded,
# {} tuple parts when the file is absent (=> exact-overlap fallback).
_graph: Optional[tuple[dict, dict, dict]] = None
# match() is gene-keyed and the pipeline calls it once per post-QC variant (~24k on an
# exome), where variants share genes heavily; memoise by (gene, patient terms). The
# loaders clear this on any (re)load so a changed graph/table never serves a stale hit.
_match_cache: dict = {}


def _read_lines(fp):
    """Yield lines from a plain or gzip-compressed TSV (the full HPO tables ship
    gzipped to keep the repo lean)."""
    if str(fp).endswith(".gz"):
        import gzip
        with gzip.open(fp, "rt") as fh:
            for line in fh:
                yield line.rstrip("\n")
    else:
        yield from fp.read_text().splitlines()


def _load() -> dict[str, set[str]]:
    global _gene_terms
    if _gene_terms is None:
        _match_cache.clear()
        d: dict[str, set[str]] = {}
        names: dict[str, str] = {}
        fp = config.HPO_GENES_LOCAL
        if fp.exists():
            for line in _read_lines(fp):
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


def _load_graph() -> tuple[dict, dict, dict]:
    global _graph
    if _graph is None:
        _match_cache.clear()
        parents: dict[str, list[str]] = {}
        ic: dict[str, float] = {}
        names: dict[str, str] = {}
        fp = config.HPO_GRAPH_LOCAL
        if fp.exists():
            for line in _read_lines(fp):
                if not line.strip() or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 3:
                    continue
                hid = parts[0]
                names[hid] = parts[1]
                try:
                    ic[hid] = float(parts[2])
                except ValueError:
                    ic[hid] = 0.0
                parents[hid] = [p for p in (parts[3].split("|") if len(parts) > 3 and parts[3] else []) if p]
        _graph = (parents, ic, names)
    return _graph


def _ancestors(term: str, parents: dict, cache: dict) -> set:
    """Transitive is_a closure incl. self. Iterative + visited-set: a DAG's multiple
    parents are fine and a corrupted cyclic graph degrades instead of overflowing."""
    if term in cache:
        return cache[term]
    seen: set = set()
    stack = [term]
    while stack:
        t = stack.pop()
        if t in seen:
            continue
        seen.add(t)
        stack.extend(p for p in parents.get(t, ()) if p not in seen)
    cache[term] = seen
    return seen


def _lin(a: str, b: str, parents: dict, ic: dict, cache: dict) -> float:
    """Lin (1998) similarity in [0,1]: 2·IC(MICA) / (IC(a)+IC(b))."""
    if a == b:
        return 1.0
    common = _ancestors(a, parents, cache) & _ancestors(b, parents, cache)
    if not common:
        return 0.0
    mica = max((ic.get(t, 0.0) for t in common), default=0.0)
    denom = ic.get(a, 0.0) + ic.get(b, 0.0)
    if denom <= 0:
        return 0.0
    return max(0.0, min(1.0, 2.0 * mica / denom))


def _semantic_match(gene_terms: set, patient_hpo: list[str]) -> Optional[dict]:
    """Ontology-aware score, or None if no graph is available (caller falls back)."""
    parents, ic, names = _load_graph()
    if not parents:
        return None
    cache: dict = {}
    per_term = []  # (patient_term, best_gene_term, sim)
    for p in patient_hpo:
        best_g, best_s = None, 0.0
        for g in gene_terms:
            s = _lin(p, g, parents, ic, cache)
            if s > best_s:
                best_s, best_g = s, g
        per_term.append((p, best_g, best_s))
    # score = best-match-average (overall phenotype coverage, drives PP4); best =
    # the single strongest patient<->gene match (drives primary-vs-secondary routing,
    # so a gene that strongly explains ONE key phenotype isn't diluted out of primary).
    score = round(sum(t[2] for t in per_term) / len(patient_hpo), 3) if patient_hpo else 0.0
    best = round(max((t[2] for t in per_term), default=0.0), 3)
    matched = [f"{p}→{g} ({names.get(g, '?')}, {s:.2f})"
               for (p, g, s) in per_term if g and s >= 0.3]
    return {"score": score, "best": best, "matched_terms": matched,
            "_source": "HPO ontology-aware (Lin/IC, local)"}


def match(gene: Optional[str], patient_hpo: list[str]) -> dict:
    """Return {'score','matched_terms','_source'} for a gene vs the patient's terms.

    ``score`` is in [0,1]. With the ontology graph it is a best-match-average Lin
    similarity (specificity-weighted, credits related terms); without it, the exact
    fraction of the patient's terms the gene is directly annotated with.
    """
    if not gene or not patient_hpo:
        return {"score": 0.0, "best": 0.0, "matched_terms": [], "_source": "HPO (no gene/terms)"}
    _load()          # loaders clear _match_cache on a (re)load, so hits are never stale
    _load_graph()
    ck = (gene, tuple(patient_hpo))
    hit = _match_cache.get(ck)
    if hit is not None:
        return hit
    r = _match_compute(gene, patient_hpo)
    _match_cache[ck] = r
    return r


def _match_compute(gene: str, patient_hpo: list[str]) -> dict:
    gene_terms = _load().get(gene, set())
    if not gene_terms:
        return {"score": 0.0, "best": 0.0, "matched_terms": [],
                "_source": "HPO genes_to_phenotype (local)"}
    sem = _semantic_match(gene_terms, patient_hpo)
    if sem is not None:
        return sem
    # Fallback: exact overlap (dependency-free). 'best' is 1.0 on any exact hit so a
    # single overlapping term still routes to primary (the pre-graph ">0" behaviour).
    patient = set(patient_hpo)
    matched = sorted(patient & gene_terms)
    score = round(len(matched) / len(patient), 3) if patient else 0.0
    labelled = [f"{t} ({_term_names.get(t, '?')})" for t in matched]
    return {"score": score, "best": 1.0 if matched else 0.0, "matched_terms": labelled,
            "_source": "HPO genes_to_phenotype (local)"}
