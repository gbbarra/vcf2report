#!/usr/bin/env bash
# Download the published, checksummed ClinVar Parquet store into the LOCAL default location
# (data/clinvar/clinvar_parquet/), where vcf2report auto-detects it — no ClinVar VCF download +
# rebuild needed. Fast alternative to build_clinvar_local.py -> build_clinvar_parquet.py.
#
#   scripts/fetch_clinvar_parquet.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$REPO_ROOT/data/clinvar"
STORE="$DEST/clinvar_parquet"
TAG="clinvar-parquet-latest"
ASSET="clinvar_parquet.tar.zst"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found — https://cli.github.com  (or build with build_clinvar_local.py -> build_clinvar_parquet.py)" >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd / apt install zstd)." >&2; exit 1; }

if [ -e "$STORE" ]; then
  echo "A store already exists at $STORE — remove it to re-fetch. Nothing to do." >&2
  exit 0
fi

echo "Downloading $TAG ($ASSET) ..." >&2
gh release download "$TAG" -D "$WORK" -p "$ASSET" -p "SHA256SUMS"
echo "Verifying checksum ..." >&2
( cd "$WORK" && shasum -a 256 -c SHA256SUMS )
echo "Extracting -> $STORE ..." >&2
mkdir -p "$DEST"
zstd -dc "$WORK/$ASSET" | tar -C "$DEST" -xf -
[ -f "$STORE/_manifest.json" ] || { echo "ERROR: extraction did not produce $STORE/_manifest.json" >&2; exit 1; }
echo "OK: $STORE is ready — vcf2report auto-detects it." >&2
python3 -c "import duckdb" 2>/dev/null || \
  echo "NOTE: the 'duckdb' package is not installed — run 'pip install duckdb' so vcf2report can read the store." >&2
