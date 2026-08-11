#!/usr/bin/env python3
"""Merge per-chromosome results.json files into one whole-genome results.json.

Used by the chunked WGS path (docs/WGS.md) after each chromosome has been
classified on its own. Concatenation of findings is exact (ACMG is per-variant);
the sequencing-quality panel is recombined with per-chunk count weighting.

    python scripts/merge_wgs_results.py OUT/chr*.results.json --out OUT/merged.results.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("inputs", nargs="+", help="per-chromosome results.json files")
    p.add_argument("--out", required=True, help="path for the merged results.json")
    args = p.parse_args()

    from vcf2report.report.merge import merge_results

    # Sort by chromosome order (natural), not the shell's alphabetical glob, so the
    # concatenated findings read chr1, chr2, ... chr10 rather than chr1, chr10, chr11.
    def _chrom_key(path: str) -> tuple:
        name = Path(path).name
        import re
        m = re.search(r"chr?(\w+)", name)
        tok = (m.group(1) if m else name).upper()
        order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
        return (order.get(tok, int(tok)) if tok.isdigit() or tok in order else 99, name)

    files = sorted(args.inputs, key=_chrom_key)
    results = [json.loads(Path(f).read_text()) for f in files]
    merged = merge_results(results)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(merged, indent=2) + "\n")
    qc = merged.get("qc", {})
    print(f"merged {len(files)} chunks -> {out}", file=sys.stderr)
    print(f"  total_variants={qc.get('total_variants')} candidates={qc.get('candidates')} "
          f"classifications={len(merged.get('classifications', []))}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
