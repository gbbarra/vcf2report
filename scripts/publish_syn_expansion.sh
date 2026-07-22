#!/usr/bin/env bash
# Publish the SYN cohort EXPANSION (SYN-101..200) as a checksummed GitHub Release asset.
#
# 100 MORE validation exomes: 100 further DISTINCT 1000G DRAGEN backgrounds (none reused from
# SYN-001..100), 100 new disease genes, variants oversampled toward the VUS-producing consequences
# (missense / in-frame) to stress the cases the engine conservatively defers on. Raw, un-annotated.
# Faithful genotypes are applied on top via scripts/build_v2_biallelic.py +
# data/synthetic_cohort/v2_faithful_plan_101_200.json.
#
#   scripts/publish_syn_expansion.sh
set -euo pipefail
export COPYFILE_DISABLE=1  # macOS: don't let tar synthesize ._* AppleDouble sidecars

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
V="$DIR/v2_build"
TAG="syn-cohort-expansion"
ASSET="syn_cohort_expansion.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }
n="$(find "$V" -name 'SYN-*.synthetic.vcf.gz' 2>/dev/null | wc -l | tr -d ' ')"
[ "$n" -ge 100 ] || { echo "ERROR: only $n/100 expansion VCFs in $V — run scripts/run_cohort_v2_loop.sh first." >&2; exit 1; }
[ -f "$V/truth.tsv" ] || { echo "ERROR: $V/truth.tsv missing." >&2; exit 1; }

echo "Materialising HPO sidecars from cohort_v2.tsv (never ship an empty --hpo file) ..." >&2
python3 "$REPO_ROOT/scripts/fill_hpo_sidecars.py" "$DIR/cohort_v2.tsv" "$V"
empty="$(find "$V" -name 'SYN-*.hpo.txt' -empty | wc -l | tr -d ' ')"
[ "$empty" -eq 0 ] || { echo "ERROR: $empty empty HPO sidecars in $V after fill — aborting publish." >&2; exit 1; }

echo "Packaging $n expansion exomes + HPO + truth -> $ASSET ..." >&2
( cd "$DIR" && tar --exclude='._*' -cf - \
    v2_build/SYN-*.synthetic.vcf.gz v2_build/SYN-*.synthetic.vcf.gz.tbi v2_build/SYN-*.hpo.txt \
    v2_build/truth.tsv cohort_v2.tsv v2_faithful_plan_101_200.json ) \
  | zstd -3 -T0 -o "$WORK/$ASSET" -q
( cd "$WORK" && shasum -a 256 "$ASSET" > SHA256SUMS )
echo "  size: $(du -h "$WORK/$ASSET" | cut -f1)" >&2
cat "$WORK/SHA256SUMS" >&2

NOTES="$WORK/notes.md"
cat > "$NOTES" <<'EOF'
100 MORE synthetic validation exomes for vcf2report (SYN-101..200) — the cohort expansion.

100 further DISTINCT 1000G DRAGEN v4.4.7 backgrounds (none reused from SYN-001..100), 100 NEW disease
genes, each spiked with one variant from a real GA4GH Phenopacket-Store case (with that case's HPO and
its faithful genotype). The variant selection deliberately OVERSAMPLES the VUS-producing consequences —
67/100 missense + in-frame (vs 45 in the first 100) — to build up N on the cases the conservative engine
defers on and to exercise the probable-pathogenic VUS triage.

Contents: v2_build/SYN-1NN.synthetic.vcf.gz(+.tbi) (raw, un-annotated), v2_build/SYN-1NN.hpo.txt,
v2_build/truth.tsv, cohort_v2.tsv (config), v2_faithful_plan_101_200.json (the faithful second alleles).

Measured (docs/BENCHMARK.md): diagnostic sensitivity 89/100 on this independent set — the same as the
first 100's 91 (~180/200 overall), confirming the number is a property of the engine, not one cohort.
See docs/COHORT_CONSTRUCTION.md for exactly how these were built. De-identified synthetic data (public
1000G backgrounds + public phenopacket variants), not real patient data, not for clinical use.

Apply faithful genotypes: scripts/build_v2_biallelic.py --plan data/synthetic_cohort/v2_faithful_plan_101_200.json \
  --cohort-tsv data/synthetic_cohort/cohort_v2.tsv --src-dir data/synthetic_cohort/v2_build --out <dir>
Fetch + verify: scripts/fetch_syn_expansion.sh
EOF

echo "Publishing release $TAG ..." >&2
if gh release view "$TAG" >/dev/null 2>&1; then
  gh release upload "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" --clobber
else
  gh release create "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" \
    --title "SYN validation cohort — expansion SYN-101..200 (vcf2report)" --notes-file "$NOTES"
fi
echo "Done. Anyone can now: scripts/fetch_syn_expansion.sh" >&2
