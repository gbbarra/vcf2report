#!/usr/bin/env python3
"""Build the compact HPO ontology graph used for ontology-aware phenotype matching.

Downloads the Human Phenotype Ontology (``hp.obo``), parses the ``is_a`` graph,
and computes each term's Information Content (IC) by propagating the bundled
gene->phenotype annotations up the graph (a gene annotated to a term is implicitly
annotated to all of that term's ancestors). Writes a small, committed file:

    data/hpo/hpo_graph.tsv.gz   columns: hpo_id  name  ic  parents(|-separated)

Runtime (``annotate/hpo.py``) loads this for a Lin/Resnik-style similarity; if the
file is absent it falls back to exact term overlap, so this build is optional.

    VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_hpo_graph.py
"""
from __future__ import annotations

import gzip
import math
import os
import sys
from pathlib import Path
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from vcf2report import config  # noqa: E402

HP_OBO_URL = "https://purl.obolibrary.org/obo/hp.obo"
ROOT = "HP:0000118"  # Phenotypic abnormality


def _require_network() -> None:
    if os.environ.get("VCF2REPORT_ALLOW_NETWORK", "").strip().lower() not in {"1", "true", "yes"}:
        sys.exit("Refusing to hit the network. Re-run with VCF2REPORT_ALLOW_NETWORK=1")


def parse_obo(text: str) -> tuple[dict[str, list[str]], dict[str, str]]:
    """Return (term -> direct parents, term -> name) from an hp.obo body."""
    parents: dict[str, list[str]] = {}
    names: dict[str, str] = {}
    cur: str | None = None
    obsolete = False
    for line in text.splitlines():
        if line == "[Term]":
            cur, obsolete = None, False
            continue
        if line.startswith("[") and line != "[Term]":
            cur = None  # a non-Term stanza (e.g. [Typedef])
            continue
        if cur is None and line.startswith("id: HP:"):
            cur = line[4:].strip()
            parents.setdefault(cur, [])
            continue
        if cur is None:
            continue
        if line.startswith("name:"):
            names[cur] = line[5:].strip()
        elif line.startswith("is_a:"):
            parents[cur].append(line[5:].split("!")[0].strip())
        elif line.startswith("is_obsolete: true"):
            obsolete = True
            parents.pop(cur, None)
            names.pop(cur, None)
            cur = None
    return parents, names


def ancestors(term: str, parents: dict[str, list[str]], cache: dict[str, set]) -> set[str]:
    if term in cache:
        return cache[term]
    anc = {term}
    for p in parents.get(term, ()):  # tolerate dangling parents
        if p in parents:
            anc |= ancestors(p, parents, cache)
        else:
            anc.add(p)
    cache[term] = anc
    return anc


def compute_ic(parents: dict[str, list[str]], gene_terms: list[str]) -> dict[str, float]:
    """IC(t) = -log( genes annotated to t-or-a-descendant / total genes )."""
    from collections import defaultdict
    ann: dict[str, set] = defaultdict(set)
    cache: dict[str, set] = {}
    # gene_terms is a flat list of (gene, hpo_id); propagate each to its ancestors.
    per_gene: dict[str, set] = defaultdict(set)
    for gene, hid in gene_terms:
        if hid in parents:
            per_gene[gene] |= ancestors(hid, parents, cache)
    total = len(per_gene) or 1
    for gene, terms in per_gene.items():
        for t in terms:
            ann[t].add(gene)
    ic: dict[str, float] = {}
    for t in parents:
        n = len(ann.get(t, ()))
        ic[t] = round(-math.log(n / total), 4) if n > 0 else 0.0
    return ic


def load_gene_terms() -> list[tuple[str, str]]:
    fp = config.HPO_GENES_LOCAL
    out: list[tuple[str, str]] = []
    with gzip.open(fp, "rt") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out.append((parts[0], parts[1]))
    return out


def main() -> None:
    _require_network()
    print(f"Downloading {HP_OBO_URL} ...", flush=True)
    with urlopen(HP_OBO_URL, timeout=120) as r:
        text = r.read().decode("utf-8", "replace")
    parents, names = parse_obo(text)
    print(f"  {len(parents)} terms parsed", flush=True)

    gene_terms = load_gene_terms()
    print(f"  {len(gene_terms)} gene->term annotations for IC", flush=True)
    ic = compute_ic(parents, gene_terms)

    out_fp = config.HPO_GRAPH_LOCAL
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_fp, "wt") as fh:
        fh.write("# HPO ontology graph for ontology-aware matching. "
                 "Columns: hpo_id\tname\tic\tparents(|)\n")
        for hid in sorted(parents):
            row = [hid, names.get(hid, ""), f"{ic.get(hid, 0.0):g}",
                   "|".join(parents.get(hid, []))]
            fh.write("\t".join(row) + "\n")
    size = out_fp.stat().st_size
    print(f"Wrote {out_fp} ({size/1024:.0f} KB, {len(parents)} terms)")


if __name__ == "__main__":
    main()
