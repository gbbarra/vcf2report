#!/usr/bin/env bash
# Auto-retry driver for the SYN cohort: S3 throughput to some networks fluctuates, so a single pass
# builds only what downloads fast and skips the slow ones. This re-runs make_syn_cohort until all
# 100 are built (resumable — skips those already done), then validates. Launch DETACHED so it
# survives closing the terminal / Claude Code, and keep the Mac awake:
#
#   cd ~/vcf2report && caffeinate -is nohup setsid bash scripts/run_cohort_loop.sh > cohort.log 2>&1 & disown
#
# Check progress:  grep -c '.' <(ls data/synthetic_cohort/SYN-*.synthetic.vcf.gz) ; tail cohort.log
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"; cd "$REPO"
DIR="data/synthetic_cohort"

count() { find "$DIR" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | wc -l | tr -d ' '; }

echo "=== cohort loop START $(date) — $(count)/100 already built ==="
while [ "$(count)" -lt 100 ]; do
  bash scripts/make_syn_cohort.sh || true
  echo "=== $(count)/100 built — $(date '+%Y-%m-%d %H:%M:%S') — S3 window pass done; retry in 3 min ==="
  [ "$(count)" -lt 100 ] && sleep 180
done
echo "=== 100/100 built — validating $(date) ==="
bash scripts/validate_cohort.sh
echo "=== COHORT COMPLETE $(date) ==="
