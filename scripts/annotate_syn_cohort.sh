#!/usr/bin/env bash
# Annotate the SYN validation cohort (gene + consequence + HGVS) — one VCF at a time.
#
# The cohort ships as raw DRAGEN calls, which carry no functional annotation. Until it is
# annotated, the engine cannot evaluate PVS1/PM4/PP3 on it, so any "no over-call" result is
# an artifact: an unannotated background CANNOT be called P/LP, no matter how the engine
# behaves. Annotating first is what makes the specificity number mean anything.
#
#   bash scripts/annotate_syn_cohort.sh            # all of SYN-001..SYN-100
#   bash scripts/annotate_syn_cohort.sh 1 10       # a range
#
# Resumable: an existing, indexed output is skipped, so an interrupted run is re-runnable.
# ~45 s and ~8.5 MB per case (~75 min / ~850 MB for 100).
set -uo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
COHORT="${VCF2REPORT_DATA:-$REPO/data}/synthetic_cohort"
DEST="$COHORT/annotated"
FROM="${1:-1}"
TO="${2:-100}"

mkdir -p "$DEST"
ok=0; skip=0; fail=0

for i in $(seq "$FROM" "$TO"); do
  id="$(printf 'SYN-%03d' "$i")"
  raw="$COHORT/$id.synthetic.vcf.gz"
  out="$DEST/$id.annotated.vcf.gz"

  [[ -f "$raw" ]] || { echo "[$id] SKIP — no raw VCF"; skip=$((skip + 1)); continue; }
  if [[ -f "$out" && -f "$out.tbi" ]]; then
    echo "[$id] skip — already annotated"; skip=$((skip + 1)); continue
  fi

  echo "[$id] annotating ..."
  if bash "$REPO/scripts/annotate_vcf.sh" "$raw" "$out" > "$DEST/$id.log" 2>&1; then
    echo "[$id] OK — $(grep -o 'annotated [0-9]*/[0-9]* records' "$DEST/$id.log" | tail -1)"
    ok=$((ok + 1))
  else
    # Never leave a half-written output behind: a truncated VCF that looks present would be
    # silently treated as annotated by the resume check on the next run.
    rm -f "$out" "$out.tbi"
    echo "[$id] FAIL — see $DEST/$id.log"; tail -3 "$DEST/$id.log" | sed 's/^/         /'
    fail=$((fail + 1))
  fi
done

echo "=== annotated: $ok · skipped: $skip · failed: $fail ==="
[[ "$fail" -eq 0 ]]
