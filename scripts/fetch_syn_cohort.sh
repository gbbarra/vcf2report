#!/usr/bin/env bash
# Download the published, checksummed SYN validation cohort (100 synthetic exomes + HPO + truth)
# into data/synthetic_cohort/ — no S3 rebuild needed.
#
#   scripts/fetch_syn_cohort.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
TAG="syn-cohort-v1"
ASSET="syn_cohort_v1.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }

if [ -n "$(find "$DIR" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | head -1)" ]; then
  echo "Cohort VCFs already present in $DIR — remove them to re-fetch. Nothing to do." >&2
  exit 0
fi

echo "Downloading $TAG ($ASSET) ..." >&2
gh release download "$TAG" -D "$WORK" -p "$ASSET" -p "SHA256SUMS"
echo "Verifying checksum ..." >&2
( cd "$WORK" && shasum -a 256 -c SHA256SUMS )
echo "Extracting -> $DIR ..." >&2
mkdir -p "$DIR"
zstd -dc "$WORK/$ASSET" | tar --exclude='._*' -C "$DIR" -xf -
echo "OK: $(find "$DIR" -name 'SYN-*.synthetic.vcf.gz' | wc -l | tr -d ' ') exomes ready in $DIR." >&2
echo "Validate: bash scripts/validate_cohort.sh" >&2
