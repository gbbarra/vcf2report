#!/usr/bin/env python3
"""Run the full pipeline headless (no Claude). Thin wrapper over the CLI.

    python scripts/run_headless.py [VCF] --hpo HPO_FILE --stdout

Defaults to the bundled sample VCF + HPO terms.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report.cli import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
