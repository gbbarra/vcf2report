#!/usr/bin/env bash
# REQUIRED install step — provision the 3 annotation Parquet stores the analysis needs, then
# verify them with the store gate. Run ONCE per machine after `git clone` + `pip install -e .`.
#
#   bash scripts/setup_stores.sh
#
# gnomAD v4.1 + ClinVar : downloaded PRE-BUILT + checksummed from the GitHub releases (fast).
# AlphaMissense         : fetched from DeepMind (CC BY-NC-SA 4.0) + built LOCALLY — it is NOT
#                         redistributed by this project; you download it under its own licence.
# Prereqs: gh (GitHub CLI), zstd, duckdb (pip install duckdb); AlphaMissense also needs htslib/tabix.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$REPO_ROOT"

echo "== [1/3] gnomAD v4.1 Parquet — pre-built release (~1.3 GB) =="
bash scripts/fetch_gnomad_parquet.sh

echo "== [2/3] ClinVar Parquet — pre-built release (~60 MB, weekly) =="
bash scripts/fetch_clinvar_parquet.sh

echo "== [3/3] AlphaMissense — CC BY-NC-SA: fetch from DeepMind + build locally (~1 GB) =="
if [ -d data/alphamissense/am_parquet ]; then
  echo "  am_parquet already present — skipping."
else
  bash scripts/fetch_alphamissense.sh
  python3 scripts/build_alphamissense_parquet.py
fi

echo ""
echo "== verifying all stores (gate) =="
python3 scripts/check_stores.py --gate
echo ""
echo "Setup complete. Try it: python3 scripts/run_headless.py   (or use the /vcf2report skill)"
