#!/usr/bin/env python3
"""Build the real HPO gene->phenotype table used for phenotype prioritisation.

The Human Phenotype Ontology distributes ``genes_to_phenotype.txt`` only as a
release artifact (GitHub release download / ontology.jax.org), neither of which
is reachable from a locked-down build network. The ``pyhpo`` PyPI package,
however, bundles the full ontology + gene annotations and installs cleanly, so we
use it as the reproducible source and export to the app's simple TSV
(``gene<TAB>hpo_id<TAB>hpo_name``) that ``annotate/hpo.py`` reads.

This is real, public HPO data — safe to commit — and replaces the tiny curated
subset previously bundled.

    pip install pyhpo
    python scripts/fetch_hpo.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "data" / "hpo" / "genes_to_phenotype.tsv"


def main() -> int:
    try:
        from pyhpo import Ontology
    except ImportError:
        print("pip install pyhpo first", file=sys.stderr)
        return 1
    Ontology()
    rows = 0
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as w:
        w.write("# HPO gene-to-phenotype (full, from pyhpo / Human Phenotype "
                "Ontology). Columns: gene\thpo_id\thpo_name\n")
        for gene in sorted(Ontology.genes, key=lambda g: g.name):
            for tid in sorted(gene.hpo):
                term = Ontology[tid]
                w.write(f"{gene.name}\t{term.id}\t{term.name}\n")
                rows += 1
    print(f"Wrote {rows} gene-phenotype rows "
          f"({len(Ontology.genes)} genes) to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
