#!/usr/bin/env python3
"""Stage 1 dependency check — prints the SAME readiness JSON as the MCP data_status
tool, so the guided flow reports identically in Claude Code (Bash) and Claude Desktop
(MCP). Emits python version, annotation tools on PATH, each local store + what it
enables, and the network-egress flag.

    python3 scripts/preflight.py
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report.status import readiness  # noqa: E402

if __name__ == "__main__":
    print(json.dumps(readiness(), indent=2))
