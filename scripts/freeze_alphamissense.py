#!/usr/bin/env python3
"""Freeze AlphaMissense scores for the concordance panel — OFFLINE (local tabix).

Reads ``data/concordance/ground_truth.tsv`` and the local AlphaMissense hg38 tabix
file (fetch once via scripts/fetch_alphamissense.sh), writing am_pathogenicity +
am_class per panel variant into ``data/concordance/alphamissense_frozen.json``.
No network — a pure local lookup, so it is safe to run any time. Idempotent /
resumable (existing entries are kept).

    python scripts/freeze_alphamissense.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import config  # noqa: E402
from vcf2report.annotate import alphamissense  # noqa: E402
from vcf2report.models import Variant  # noqa: E402


def main() -> int:
    gt = config.DATA_DIR / "concordance" / "ground_truth.tsv"
    out = config.DATA_DIR / "concordance" / "alphamissense_frozen.json"

    if not config.ALPHAMISSENSE_LOCAL.exists():
        print(f"ERROR: {config.ALPHAMISSENSE_LOCAL} not found.\n"
              f"  Run: bash scripts/fetch_alphamissense.sh", file=sys.stderr)
        return 2
    if not gt.exists():
        print("ERROR: ground_truth.tsv not found; build the panel first "
              "(scripts/build_concordance_panel.py).", file=sys.stderr)
        return 2

    keys = [line.split("\t")[0] for line in gt.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    frozen: dict = json.loads(out.read_text()) if out.exists() else {}

    todo = [k for k in keys if k not in frozen]
    print(f"Freezing AlphaMissense for {len(todo)} variants "
          f"({len(frozen)} already cached)...")
    for i, key in enumerate(todo, 1):
        chrom, pos, ref, alt = key.split("-")
        r = alphamissense.lookup(Variant(chrom=chrom, pos=int(pos), ref=ref, alt=alt))
        frozen[key] = {"am_pathogenicity": r.get("am_pathogenicity"),
                       "am_class": r.get("am_class")}
        if i % 25 == 0 or i == len(todo):
            out.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")

    out.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")
    scored = sum(1 for v in frozen.values() if v.get("am_pathogenicity") is not None)
    print(f"Done: {len(frozen)} variants frozen, {scored} with a missense "
          f"AlphaMissense score (non-missense variants have none).")
    print(f"  {out}")
    print("Now run:  python scripts/run_concordance.py   (offline)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
