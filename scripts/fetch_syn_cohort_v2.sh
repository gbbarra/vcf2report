#!/usr/bin/env bash
# Download the published, checksummed FAITHFUL biallelic SYN cohort (v2) into
# data/synthetic_cohort/v2/ (+ shared HPO/truth/cohort in data/synthetic_cohort/). No rebuild needed.
#
#   scripts/fetch_syn_cohort_v2.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
TAG="syn-cohort-v2"
ASSET="syn_cohort_v2.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }

if [ -n "$(find "$DIR/v2" -name 'SYN-*.v2.vcf.gz' 2>/dev/null | head -1)" ]; then
  echo "v2 VCFs already present in $DIR/v2 — remove them to re-fetch. Nothing to do." >&2
  exit 0
fi

echo "Downloading $TAG ($ASSET) ..." >&2
gh release download "$TAG" -D "$WORK" -p "$ASSET" -p "SHA256SUMS"
echo "Verifying checksum ..." >&2
( cd "$WORK" && shasum -a 256 -c SHA256SUMS )
echo "Extracting -> $DIR ..." >&2
mkdir -p "$DIR"
zstd -dc "$WORK/$ASSET" | tar -C "$DIR" -xf -
echo "OK: $(find "$DIR/v2" -name 'SYN-*.v2.vcf.gz' | wc -l | tr -d ' ') faithful exomes ready in $DIR/v2." >&2
echo "Annotate: bash scripts/annotate_syn_cohort.sh   (point it at the v2 dir)" >&2
