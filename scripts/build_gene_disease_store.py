#!/usr/bin/env python3
"""Freeze the gene -> diseases (+ per-disease inheritance) store for the laudo.

Joins the HPO project's two release files into ``data/hpo/gene_disease.json.gz``
(read at runtime by ``vcf2report.annotate.gene_disease``):

  * ``genes_to_phenotype.txt`` — gene x disease x HPO term. Inheritance is itself an
    HPO term (HP:0000006 AD, HP:0000007 AR, ...), so grouping by (gene, disease_id) and
    keeping the inheritance terms gives inheritance **per disease**, not per gene.
  * ``phenotype.hpoa`` — disease_id -> readable disease_name.

Both from https://github.com/obophenotype/human-phenotype-ontology/releases (same release).

    python scripts/build_gene_disease_store.py genes_to_phenotype.txt phenotype.hpoa \
        --out data/hpo/gene_disease.json.gz
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def _header_index(line: str, wanted: list[str]) -> dict[str, int]:
    cols = line.lstrip("#").rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(cols)}
    return {w: idx[w] for w in wanted if w in idx}


def _disease_names(hpoa: Path) -> dict[str, str]:
    """phenotype.hpoa -> {database_id: disease_name} (first seen)."""
    names: dict[str, str] = {}
    with open(hpoa) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            if header is None:
                header = _header_index(line, ["database_id", "disease_name"])
                continue
            f = line.rstrip("\n").split("\t")
            try:
                did, name = f[header["database_id"]], f[header["disease_name"]]
            except (KeyError, IndexError):
                continue
            names.setdefault(did, name)
    return names


def main() -> int:
    from vcf2report.annotate.gene_disease import INHERITANCE_TERMS

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genes_to_phenotype")
    p.add_argument("phenotype_hpoa")
    p.add_argument("--out", default="data/hpo/gene_disease.json.gz")
    args = p.parse_args()

    names = _disease_names(Path(args.phenotype_hpoa))
    print(f"loaded {len(names)} disease names from {args.phenotype_hpoa}", file=sys.stderr)

    # (gene, disease_id) -> set of inheritance labels
    inh: dict[tuple, set] = {}
    with open(args.genes_to_phenotype) as fh:
        header = None
        for line in fh:
            if header is None:
                header = _header_index(line, ["gene_symbol", "hpo_id", "disease_id"])
                if not {"gene_symbol", "hpo_id", "disease_id"} <= header.keys():
                    raise SystemExit(f"unexpected columns in {args.genes_to_phenotype}: "
                                     "need gene_symbol, hpo_id, disease_id")
                continue
            f = line.rstrip("\n").split("\t")
            try:
                gene = f[header["gene_symbol"]].upper()
                hpo_id = f[header["hpo_id"]]
                did = f[header["disease_id"]]
            except IndexError:
                continue
            if not gene or not did:
                continue
            key = (gene, did)
            inh.setdefault(key, set())
            if hpo_id in INHERITANCE_TERMS:
                inh[key].add(INHERITANCE_TERMS[hpo_id])

    store: dict[str, list] = {}
    for (gene, did), labels in inh.items():
        store.setdefault(gene, []).append({
            "disease_id": did,
            "disease_name": names.get(did, did),
            "inheritance": sorted(labels),
        })
    for gene in store:
        store[gene].sort(key=lambda d: d["disease_name"])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt") as fh:
        json.dump(store, fh)
    n_dis = sum(len(v) for v in store.values())
    print(f"wrote {out}: {len(store)} genes, {n_dis} gene-disease links", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
