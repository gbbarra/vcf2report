#!/usr/bin/env bash
# Refresh the ClinGen Haploinsufficiency=3 gene list (the PVS1 curated LoF-mechanism route).
# Public ClinGen Dosage Sensitivity data. Committed to the repo (tiny); re-run to update.
#   bash scripts/fetch_clingen_hi.sh
set -euo pipefail
REPO="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$REPO/data/constraint/clingen_haploinsufficiency.tsv"
URL="https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv"
mkdir -p "$(dirname "$OUT")"
tmp="$(mktemp)"; trap 'rm -f "$tmp"' EXIT
curl -fsSL --max-time 120 -o "$tmp" "$URL"
{ echo "# ClinGen Dosage Sensitivity — genes with Haploinsufficiency=3 (sufficient evidence that"
  echo "# loss of function causes disease). Source: ftp.clinicalgenome.org (public). One gene per line."
  awk -F'\t' '!/^#/ && $1!="Gene Symbol" && $5==3 {print $1}' "$tmp" | sort -u
} > "$OUT"
echo "Wrote $(grep -vc '^#' "$OUT") HI=3 genes -> $OUT"
