#!/usr/bin/env bash
# Publish the locally-built gnomAD Parquet store as a checksummed GitHub Release asset,
# so the online jury (and anyone) can fetch + verify the exact store without the ~150 GB
# from-scratch build. Run once, after scripts/build_gnomad_parquet.py.
#
#   scripts/publish_gnomad_parquet.sh [STORE_DIR]   # default: data/gnomad/gnomad_parquet
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STORE="${1:-$REPO_ROOT/data/gnomad/gnomad_parquet}"
TAG="gnomad-parquet-v4.1"
ASSET="gnomad_parquet_v4.1.tar.zst"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

[ -d "$STORE" ] || { echo "ERROR: store not found: $STORE — build it first with scripts/build_gnomad_parquet.py" >&2; exit 1; }
[ -f "$STORE/_meta.json" ] || { echo "ERROR: $STORE/_meta.json missing — the build looks incomplete." >&2; exit 1; }
command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found — https://cli.github.com" >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd / apt install zstd)." >&2; exit 1; }

echo "Packaging $STORE -> $ASSET ..." >&2
tar -C "$(dirname "$STORE")" -cf - "$(basename "$STORE")" | zstd -19 -T0 -o "$WORK/$ASSET" -q
( cd "$WORK" && shasum -a 256 "$ASSET" > SHA256SUMS )
echo "  size: $(du -h "$WORK/$ASSET" | cut -f1)" >&2
cat "$WORK/SHA256SUMS" >&2

echo "Publishing release $TAG ..." >&2
if gh release view "$TAG" >/dev/null 2>&1; then
  gh release upload "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" --clobber
else
  gh release create "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" \
    --title "gnomAD v4.1 frequency Parquet (vcf2report)" \
    --notes-file "$REPO_ROOT/data/gnomad/NOTICE.md"
fi
echo "Done. Anyone can now: scripts/fetch_gnomad_parquet.sh" >&2
