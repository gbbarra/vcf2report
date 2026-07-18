#!/usr/bin/env bash
# Detached, resumable builder for the SYN cohort v2 EXPANSION (SYN-101..200): 100 NEW distinct
# 1000G backgrounds + phenopacket variants oversampled toward missense/inframe (the VUS-producing
# consequences). Streams DRAGEN VCFs from S3 (throttled — hence the retry loop) and skips any case
# already built, so it survives interruptions. Run detached; it self-retries until 100 are built.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"; cd "$REPO"
COHORT="$REPO/data/synthetic_cohort/cohort_v2.tsv"
OUT="$REPO/data/synthetic_cohort/v2_build"
count() { find "$OUT" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | wc -l | tr -d ' '; }
echo "=== v2 cohort loop START $(date) — $(count)/100 already built ==="
tries=0
while [ "$(count)" -lt 100 ] && [ "$tries" -lt 200 ]; do
  bash scripts/make_syn_cohort.sh "$COHORT" "$OUT" || true
  tries=$((tries+1))
  echo "=== $(count)/100 built — $(date '+%H:%M:%S') — pass $tries done; retry in 3 min ==="
  [ "$(count)" -lt 100 ] && sleep 180
done
echo "=== v2 build: $(count)/100 built — $(date) ==="
