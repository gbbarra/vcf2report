#!/usr/bin/env bash
# Publish the FAITHFUL biallelic SYN cohort (v2) as a checksummed GitHub Release asset.
#
# v2 differs from v1 in one way: each planted variant was matched back to its GA4GH
# Phenopacket-Store case and the patient's REAL genotype restored — 22 compound-heterozygous
# (both true alleles), 35 homozygous (the real allelicState), 43 single-allele (the source
# genuinely recorded one). v1 planted every variant as a lone heterozygote, which silently
# turned every recessive patient into a carrier and could not test recessive diagnostic recovery.
# Recipe (committed, reproducible): scripts/build_v2_biallelic.py + data/synthetic_cohort/v2_faithful_plan.json.
#
#   scripts/publish_syn_cohort_v2.sh
set -euo pipefail
export COPYFILE_DISABLE=1  # macOS: don't let tar synthesize ._* AppleDouble sidecars

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$REPO_ROOT/data/synthetic_cohort"
V2="$DIR/v2"
TAG="syn-cohort-v2"
ASSET="syn_cohort_v2.tar.zst"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found." >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd)." >&2; exit 1; }
n="$(find "$V2" -name 'SYN-*.v2.vcf.gz' 2>/dev/null | wc -l | tr -d ' ')"
[ "$n" -ge 100 ] || { echo "ERROR: only $n/100 v2 VCFs in $V2 — run scripts/build_v2_biallelic.py --out $V2 first." >&2; exit 1; }
[ -f "$V2/v2_faithful_plan.json" ] || { echo "ERROR: $V2/v2_faithful_plan.json missing." >&2; exit 1; }

echo "Materialising HPO sidecars from cohort.tsv, co-located in v2/ (never ship an empty --hpo file) ..." >&2
python3 "$REPO_ROOT/scripts/fill_hpo_sidecars.py" "$DIR/cohort.tsv" "$V2"
empty="$(find "$V2" -name 'SYN-*.hpo.txt' -empty | wc -l | tr -d ' ')"
[ "$empty" -eq 0 ] || { echo "ERROR: $empty empty HPO sidecars in $V2 after fill — aborting publish." >&2; exit 1; }

echo "Packaging $n faithful biallelic exomes + HPO + truth -> $ASSET ..." >&2
( cd "$DIR" && tar --exclude='._*' -cf - \
    v2/SYN-*.v2.vcf.gz v2/SYN-*.v2.vcf.gz.tbi v2/SYN-*.hpo.txt \
    v2/v2_manifest.json v2/v2_faithful_plan.json truth.tsv cohort.tsv ) \
  | zstd -3 -T0 -o "$WORK/$ASSET" -q
( cd "$WORK" && shasum -a 256 "$ASSET" > SHA256SUMS )
echo "  size: $(du -h "$WORK/$ASSET" | cut -f1)" >&2
cat "$WORK/SHA256SUMS" >&2

NOTES="$WORK/notes.md"
cat > "$NOTES" <<'EOF'
100 synthetic validation exomes for vcf2report — **v2, faithful to the source phenopacket genotype**.

Same 100 distinct 1000G DRAGEN v4.4.7 backgrounds (MANE/GENCODE exome BED, ~100k variants) and the
same 100 distinct causative genes as v1, but each planted variant now carries the patient's REAL
genotype, recovered by matching to its GA4GH Phenopacket-Store case:
  * 22 compound-heterozygous (both true alleles)
  * 35 homozygous (the real allelicState)
  * 43 single-allele (the source genuinely recorded one)

Why: v1 planted every variant as a lone heterozygote, which for an autosomal-recessive gene is a
healthy CARRIER, not a diagnosis — so v1 could not test recessive diagnostic recovery. v2 fixes that.

Contents: v2/SYN-00N.v2.vcf.gz(+.tbi) (raw, un-annotated — annotate with scripts/annotate_syn_cohort.sh),
v2/SYN-00N.hpo.txt (the case's HPO terms, one per line, co-located with each VCF),
v2/v2_manifest.json (per-case genotype mode), v2/v2_faithful_plan.json (the second alleles),
truth.tsv, cohort.tsv.

Measured recovery (docs/BENCHMARK.md): diagnostic sensitivity 91/100 (compound-het 22/22, hom 32/35);
the 9 non-primary are honest limitations (non-coding RNA, HPO-dropped genes, missense VUS, one
sub-threshold phenotype), not engine misses. De-identified synthetic data (public 1000G backgrounds +
public phenopacket variants), not real patient data, not for clinical use.

Fetch + verify: scripts/fetch_syn_cohort_v2.sh
EOF

echo "Publishing release $TAG ..." >&2
if gh release view "$TAG" >/dev/null 2>&1; then
  gh release upload "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" --clobber
else
  gh release create "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" \
    --title "SYN validation cohort v2 — faithful biallelic (vcf2report)" --notes-file "$NOTES"
fi
echo "Done. Anyone can now: scripts/fetch_syn_cohort_v2.sh" >&2
