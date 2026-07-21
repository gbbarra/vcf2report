#!/usr/bin/env python3
"""Materialise per-case HPO sidecar files from a cohort TSV.

The engine takes phenotype as a separate ``--hpo`` file (one ``HP:xxxxxxx`` per
line); the VCF never carries HPO. Cohort builders record the authoritative terms
in the TSV's ``hpo`` column, but the expansion build left the ``SYN-NNN.hpo.txt``
sidecars as empty stubs — an empty file passed to ``--hpo`` yields a silent
genotype-only run (no PP4). This writes the sidecars from the TSV so no empty
HPO file is left lying next to a VCF.

    python3 scripts/fill_hpo_sidecars.py <cohort.tsv> <out_dir>

Idempotent. Reports any row whose HPO column is empty (left untouched, flagged).
"""
import csv
import sys
from pathlib import Path


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(__doc__)
        return 2
    tsv, out_dir = Path(argv[1]), Path(argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    written = empty = 0
    with tsv.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row.get("syn_id") or row.get("id")
            if not sid:
                continue
            terms = [t.strip() for t in (row.get("hpo") or "").split(",") if t.strip()]
            dest = out_dir / f"{sid}.hpo.txt"
            if not terms:
                empty += 1
                print(f"  WARN {sid}: no HPO in TSV — left {dest.name} untouched")
                continue
            dest.write_text("\n".join(terms) + "\n")
            written += 1
    print(f"wrote {written} sidecars into {out_dir} ({empty} rows had no HPO)")
    return 1 if empty else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
