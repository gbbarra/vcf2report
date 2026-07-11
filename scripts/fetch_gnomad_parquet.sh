#!/usr/bin/env bash
# Download the published, checksummed gnomAD Parquet store into the LOCAL default
# location (data/gnomad/gnomad_parquet/), where vcf2report auto-detects it — no ~150 GB
# from-scratch build needed. This is the fast alternative to build_gnomad_parquet.py.
#
#   scripts/fetch_gnomad_parquet.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="$REPO_ROOT/data/gnomad"
STORE="$DEST/gnomad_parquet"
TAG="gnomad-parquet-v4.1"
ASSET="gnomad_parquet_v4.1.tar.zst"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found — https://cli.github.com  (or build with scripts/build_gnomad_parquet.py)" >&2; exit 1; }
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
[ -f "$STORE/_meta.json" ] || { echo "ERROR: extraction did not produce $STORE/_meta.json" >&2; exit 1; }
echo "OK: $STORE is ready — vcf2report auto-detects it (no env var needed)." >&2
# Reading the store needs duckdb; fetch itself does not. Warn (don't fail) if it is absent.
python3 -c "import duckdb" 2>/dev/null || \
  echo "NOTE: the 'duckdb' package is not installed — run 'pip install duckdb' (or 'pip install .[parquet]') so vcf2report can read the store." >&2
