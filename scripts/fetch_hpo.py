#!/usr/bin/env python3
"""Build the real HPO gene->phenotype table used for phenotype prioritisation.

Preferred source is the official ``genes_to_phenotype.txt`` release artifact
(github.com/obophenotype/human-phenotype-ontology), which tracks the current HPO
annotation. The bundled ``pyhpo`` PyPI package is the offline fallback, but it ships a
frozen (older) ontology snapshot — measurably behind: the release covered 36 disease
genes the pyhpo build missed, all of them causative genes in the validation cohort, so
the stale table silently dropped their phenotype match AND their inheritance (which the
PVS1 recessive-LoF route and the PM2/BS1 ceilings read from the same table).

Exports the app's simple TSV (``gene<TAB>hpo_id<TAB>hpo_name``) that ``annotate/hpo.py``
reads, de-duplicated across diseases. HPO data is public and redistributable (HPO license).

    VCF2REPORT_ALLOW_NETWORK=1 python scripts/fetch_hpo.py        # latest release
    python scripts/fetch_hpo.py --pyhpo                           # offline fallback

The ontology graph (``hpo_graph.tsv.gz``, Lin/IC similarity) is a superset of the term
space and does not need rebuilding when this table is refreshed; verify with
``scripts/build_hpo_graph.py`` if a future release adds brand-new term ids.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "hpo" / "genes_to_phenotype.tsv.gz"
RELEASE_URL = ("https://github.com/obophenotype/human-phenotype-ontology/"
               "releases/latest/download/genes_to_phenotype.txt")
_HEADER = ("# HPO gene-to-phenotype (from the HPO genes_to_phenotype.txt release). "
           "Columns: gene\thpo_id\thpo_name\n")


def _write(pairs, genes: int) -> int:
    """pairs: iterable of (gene, hpo_id, hpo_name), de-duplicated by the caller."""
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = 0
    with gzip.open(OUT, "wt") as w:
        w.write(_HEADER)
        for g, tid, name in pairs:
            w.write(f"{g}\t{tid}\t{name}\n")
            rows += 1
    print(f"Wrote {rows} gene-phenotype rows ({genes} genes) to {OUT}")
    return rows


def from_release() -> int:
    with urllib.request.urlopen(RELEASE_URL, timeout=180) as r:
        text = r.read().decode("utf-8")
    seen: set = set()
    pairs = []
    genes: set = set()
    for line in text.splitlines():
        if not line or line.startswith("ncbi_gene_id") or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        gene, tid, name = f[1], f[2], f[3]
        genes.add(gene)
        key = (gene, tid)
        if key in seen:               # same gene-term repeats across diseases (disease_id col)
            continue
        seen.add(key)
        pairs.append((gene, tid, name))
    pairs.sort(key=lambda t: (t[0], t[1]))
    return _write(pairs, len(genes))


def from_pyhpo() -> int:
    try:
        from pyhpo import Ontology
    except ImportError:
        print("pip install pyhpo (or run with network for the release)", file=sys.stderr)
        return 1
    Ontology()
    pairs = []
    for gene in sorted(Ontology.genes, key=lambda g: g.name):
        for tid in sorted(gene.hpo):
            term = Ontology[tid]
            pairs.append((gene.name, term.id, term.name))
    return 0 if _write(pairs, len(Ontology.genes)) else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pyhpo", action="store_true", help="use the offline pyhpo snapshot")
    args = ap.parse_args()
    if args.pyhpo:
        return 0 if from_pyhpo() else 1
    try:
        from_release()
        return 0
    except Exception as exc:              # network blocked / release moved → offline fallback
        print(f"release download failed ({exc}); falling back to pyhpo", file=sys.stderr)
        return 0 if from_pyhpo() else 1


if __name__ == "__main__":
    raise SystemExit(main())
