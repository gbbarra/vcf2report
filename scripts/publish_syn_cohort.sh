#!/usr/bin/env bash
# Publish the built SYN validation cohort (100 synthetic exomes + HPO + ground truth) as a
# checksummed GitHub Release asset, so anyone can fetch the exact corpus without rebuilding it from
# S3. Run after the cohort is built (scripts/make_syn_cohort.sh).
#
#   scripts/publish_syn_cohort.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
TAG="syn-cohort-v1"
ASSET="syn_cohort_v1.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }
n="$(find "$DIR" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | wc -l | tr -d ' ')"
[ "$n" -ge 100 ] || { echo "ERROR: only $n/100 cohort VCFs built — run scripts/make_syn_cohort.sh first." >&2; exit 1; }
[ -f "$DIR/truth.tsv" ] || { echo "ERROR: $DIR/truth.tsv missing." >&2; exit 1; }

echo "Packaging $n exomes + HPO + truth -> $ASSET ..." >&2
( cd "$DIR" && tar --exclude='._*' -cf - \
    SYN-*.synthetic.vcf.gz SYN-*.synthetic.vcf.gz.tbi SYN-*.hpo.txt truth.tsv cohort.tsv ) \
  | zstd -3 -T0 -o "$WORK/$ASSET" -q
( cd "$WORK" && shasum -a 256 "$ASSET" > SHA256SUMS )
echo "  size: $(du -h "$WORK/$ASSET" | cut -f1)" >&2
cat "$WORK/SHA256SUMS" >&2

NOTES="$WORK/notes.md"
cat > "$NOTES" <<'EOF'
100 synthetic validation exomes for vcf2report: 100 DISTINCT 1000G DRAGEN v4.4.7 samples, each subset
to the vendor-neutral MANE/GENCODE exome BED (~100k variants) and spiked with ONE distinct pathogenic
variant from a real GA4GH phenopacket case (carrying that case's HPO). No repeated background, no
repeated variant. SYN-001 = NA12878 (GIAB).

Contents: SYN-00N.synthetic.vcf.gz(+.tbi), SYN-00N.hpo.txt, truth.tsv (ground truth), cohort.tsv (config).

Validation (docs): all 100 structurally correct; precise recovery 59/100 surfaced (51 engine-only,
anti-circular); over-call median 1 P/LP per case; decoy phenotype 15%. De-identified, not real patient
data, not for clinical use.

Fetch + verify: scripts/fetch_syn_cohort.sh
EOF

echo "Publishing release $TAG ..." >&2
if gh release view "$TAG" >/dev/null 2>&1; then
  gh release upload "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" --clobber
else
  gh release create "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" \
    --title "SYN validation cohort v1 (vcf2report)" --notes-file "$NOTES"
fi
echo "Done. Anyone can now: scripts/fetch_syn_cohort.sh" >&2
