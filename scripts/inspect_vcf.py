#!/usr/bin/env python3
"""Stages 3 & 5 for the Bash path — detect annotation + the ACMG capability map for a
VCF, the same data the MCP inspect_vcf / analysis_capabilities tools return.

    python3 scripts/inspect_vcf.py VCF [--hpo]

--hpo signals that phenotype terms were supplied (affects the PP4 capability).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report.inspect import analysis_capabilities, inspect_vcf  # noqa: E402

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(json.dumps({"error": "usage: inspect_vcf.py VCF [--hpo]"}))
        sys.exit(2)
    vcf = sys.argv[1]
    hpo_given = "--hpo" in sys.argv[2:]
    insp = inspect_vcf(vcf)
    caps = analysis_capabilities(vcf, hpo_given=hpo_given, inspection=insp)
    print(json.dumps({"inspect": insp, "capabilities": caps}, indent=2))
