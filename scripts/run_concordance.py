#!/usr/bin/env python3
"""Print the ClinVar-vs-engine concordance panel — fully offline.

Reads the frozen panel (``data/concordance/ground_truth.tsv`` +
``gnomad_frozen.json``, produced once by ``build_concordance_panel.py``), runs the
ACMG engine on every variant with ClinVar withheld, and prints the confusion
matrix + metrics as Markdown.

    python scripts/run_concordance.py [--json] [--out data/out/concordance.md]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import concordance, config  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description="Render the concordance panel (offline).")
    ap.add_argument("--json", action="store_true", help="emit JSON instead of Markdown")
    ap.add_argument("--out", default="", help="also write the Markdown to this path")
    ap.add_argument("--include-clinvar", action="store_true",
                    help="do NOT withhold ClinVar (measures full production behaviour, circular)")
    args = ap.parse_args()

    if not concordance.GROUND_TRUTH.exists():
        print("ERROR: panel not built yet. Run:\n"
              "  VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_concordance_panel.py",
              file=sys.stderr)
        return 2

    entries = concordance.load_panel()
    if not entries:
        print("ERROR: ground_truth.tsv has no PATH/BEN-labelled variants.", file=sys.stderr)
        return 2

    result = concordance.evaluate_panel(entries, withhold_clinvar=not args.include_clinvar)

    if args.json:
        print(json.dumps(result.to_dict(), indent=2))
    else:
        print(result.to_markdown())

    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(result.to_markdown() + "\n")
        print(f"\n(written to {out})", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
