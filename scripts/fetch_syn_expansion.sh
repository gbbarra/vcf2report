#!/usr/bin/env bash
# Download the published, checksummed SYN cohort expansion (SYN-101..200) into
# data/synthetic_cohort/v2_build/ (+ cohort_v2.tsv + faithful plan). No S3 rebuild needed.
#
#   scripts/fetch_syn_expansion.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
TAG="syn-cohort-expansion"
ASSET="syn_cohort_expansion.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }

if [ -n "$(find "$DIR/v2_build" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | head -1)" ]; then
  echo "Expansion VCFs already present in $DIR/v2_build — remove them to re-fetch. Nothing to do." >&2
  exit 0
fi

echo "Downloading $TAG ($ASSET) ..." >&2
gh release download "$TAG" -D "$WORK" -p "$ASSET" -p "SHA256SUMS"
echo "Verifying checksum ..." >&2
( cd "$WORK" && shasum -a 256 -c SHA256SUMS )
echo "Extracting -> $DIR ..." >&2
mkdir -p "$DIR"
zstd -dc "$WORK/$ASSET" | tar -C "$DIR" -xf -
echo "OK: $(find "$DIR/v2_build" -name 'SYN-*.synthetic.vcf.gz' | wc -l | tr -d ' ') expansion exomes ready in $DIR/v2_build." >&2
