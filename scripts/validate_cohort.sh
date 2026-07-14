#!/usr/bin/env bash
# Validate the SYN cohort: run each SYN-00N through the engine and check whether the PLANTED gene
# (from truth.tsv) is surfaced — engine P/LP tier, a phenotype-matched "likely explanatory" finding,
# or a >=2-star ClinVar flag. Prints a per-case result + the overall surfaced rate.
#
#   bash scripts/validate_cohort.sh [COHORT_DIR]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
DIR="${1:-$REPO/data/synthetic_cohort}"
TRUTH="$DIR/truth.tsv"
[ -f "$TRUTH" ] || { echo "ERROR: $TRUTH not found — run make_syn_cohort.sh first" >&2; exit 1; }

RES="$DIR/validation_results.tsv"; echo -e "syn_id\tgene\tsurfaced" > "$RES"
tail -n +2 "$TRUTH" | while IFS=$'\t' read -r syn sample gene variant cons hn vc disease; do
  vcf="$DIR/$syn.synthetic.vcf.gz"; hpo="$DIR/$syn.hpo.txt"
  [ -f "$vcf" ] || continue
  out="$(python3 "$REPO/scripts/run_headless.py" "$vcf" --hpo "$hpo" --sample-id "$syn" --stdout 2>/dev/null || true)"
  if echo "$out" | grep -qiE "(likely explanatory|Classified Pathogenic).*\b$gene\b|\b$gene\b.*(Pathogenic|Likely Pathogenic)"; then
    echo "  [OK ] $syn  $gene  surfaced"; echo -e "$syn\t$gene\tyes" >> "$RES"
  else
    echo "  [ - ] $syn  $gene  not surfaced"; echo -e "$syn\t$gene\tno" >> "$RES"
  fi
done

python3 - "$RES" <<'PY'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
s = sum(r["surfaced"] == "yes" for r in rows)
print(f"\nSurfaced {s}/{len(rows)} planted variants ({100*s/max(1,len(rows)):.0f}%). Results: {sys.argv[1]}")
PY
