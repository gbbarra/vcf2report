#!/usr/bin/env bash
# Fetch the hpo-spiked-exomes benchmark (200 tell-free spiked exomes) and print the scoring command.
#
# Clones (or updates) the benchmark repo, downloads its release assets via the repo's OWN fetch.sh
# (raw + SnpEff-annotated VCFs), then shows how to score vcf2report against it with run_benchmark.py.
#
#   bash scripts/fetch_benchmark.sh [DEST]     # DEST defaults to ~/hpo-spiked-exomes
#
# Requires: git, zstd, ~3 GB disk. Score where the vcf2report data stores live (check_stores.py green)
# — the benchmark VCFs are tell-free (no baked gnomAD/ClinVar INFO), so the engine must look them up.
set -euo pipefail
DEST="${1:-$HOME/hpo-spiked-exomes}"
REPO="https://github.com/gbbarra/hpo-spiked-exomes.git"

if [ -d "$DEST/.git" ]; then
  echo "Updating $DEST ..." >&2
  git -C "$DEST" pull --ff-only
else
  git clone "$REPO" "$DEST"
fi

bash "$DEST/fetch.sh"

cat >&2 <<EOF

Benchmark ready in $DEST. Score vcf2report against it (run where the data stores are present):

  python3 scripts/run_benchmark.py \\
      --annotated $DEST/realistic_annotated \\
      --bench     $DEST \\
      --out       benchmark_results.tsv --jobs 4
EOF
