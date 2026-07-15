#!/usr/bin/env bash
# Validate the SYN cohort PRECISELY: run the engine on each case and check whether the PLANTED gene
# is surfaced — engine P/LP tier, phenotype-matched primary (P/LP), or a >=2-star ClinVar flag —
# using the classification objects (NOT a markdown grep, which over-counts). Also prints the
# ENGINE-ONLY rate (excludes the ClinVar flag = the anti-circular metric).
#
#   bash scripts/validate_cohort.sh [COHORT_DIR]
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"; cd "$REPO"
DIR="${1:-data/synthetic_cohort}"
python3 scripts/adversarial_cohort.py 2>/dev/null | grep -A20 "ADVERSARIAL SUMMARY" || \
  python3 scripts/adversarial_cohort.py
