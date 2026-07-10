#!/usr/bin/env bash
# Fetch + tabix-index the AlphaMissense hg38 missense predictions (CC BY 4.0).
#
# One-time download (~1 GB) used by the AlphaMissense annotation client and the
# concordance-panel freeze. AlphaMissense predictions are © Google DeepMind,
# licensed CC BY 4.0 (attribution required) — see docs/DISCLAIMERS.md.
#
#   bash scripts/fetch_alphamissense.sh
#
# Requires: wget (or curl), htslib's bgzip + tabix (brew install htslib).
set -euo pipefail

DEST="${VCF2REPORT_DATA:-$(cd "$(dirname "$0")/.." && pwd)/data}/alphamissense"
URL="https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"
FILE="$DEST/AlphaMissense_hg38.tsv.gz"

mkdir -p "$DEST"

echo "==> Downloading AlphaMissense hg38 (~1 GB, resumable) to $FILE"
if command -v wget >/dev/null 2>&1; then
  wget -c -O "$FILE" "$URL"
else
  curl -L -C - -o "$FILE" "$URL"
fi

command -v tabix >/dev/null 2>&1 || { echo "ERROR: tabix not found (brew install htslib)"; exit 1; }

echo "==> Indexing for per-variant lookup"
# Columns: #CHROM POS REF ALT genome uniprot_id transcript_id protein_variant
#          am_pathogenicity am_class. Position col 2; '#'-prefixed lines are
#          comments (incl. the #CHROM header) and skipped by tabix automatically.
if ! tabix -s 1 -b 2 -e 2 -f "$FILE" 2>/dev/null; then
  echo "   direct index failed (file not bgzip-sorted); re-compressing with bgzip"
  zcat "$FILE" | bgzip > "$FILE.tmp" && mv "$FILE.tmp" "$FILE"
  tabix -s 1 -b 2 -e 2 -f "$FILE"
fi

echo "==> Done: $FILE (+ .tbi)"
echo "    Now the AlphaMissense client + concordance freeze can read it offline."
