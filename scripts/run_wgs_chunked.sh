#!/usr/bin/env bash
# Whole-genome (WGS) analysis in a small RAM budget — the chunked path.
#
# Running the standard pipeline on a ~5M-variant WGS VCF OOMs an 8 GB machine in two
# places: SnpEff's fixed heap, and parse_vcf loading every variant into one Python list.
# This driver sidesteps both by processing ONE chromosome at a time — annotate → classify
# → per-chrom results.json — then stitches the results back together. ACMG is per-variant,
# so per-chromosome classification is exact; only the seq-quality panel needs count-weighted
# recombination (scripts/merge_wgs_results.py). Peak RSS measured ~1.5 GB (annotate) and
# ~515 MB (classify) on the largest chromosome.
#
# Usage:
#   scripts/run_wgs_chunked.sh RAW.vcf.gz REF_GRCh38.fa HPO.txt OUTDIR
#
# Env: SNPEFF_XMX (default 3g — WGS-safe), CHROMS (override the chromosome list).
set -euo pipefail

RAW="${1:?raw WGS VCF (bgzipped) required}"
REF="${2:?GRCh38 reference FASTA required}"
HPO="${3:?patient HPO term file required}"
OUT="${4:?output directory required}"
HERE="$(cd "$(dirname "$0")" && pwd)"
export SNPEFF_XMX="${SNPEFF_XMX:-3g}"   # WGS-safe heap for scripts/annotate_vcf.sh

for tool in bcftools tabix python; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH" >&2; exit 1; }
done
mkdir -p "$OUT"
[[ -f "$RAW.tbi" || -f "$RAW.csi" ]] || { echo "[index] tabix $RAW"; tabix -p vcf "$RAW"; }

# The 25 primary contigs cover ~98.6% of a typical WGS VCF's variants; the ~195 alt/decoy/
# unplaced/HLA contigs are not worth the per-chunk overhead. Match the VCF's naming (chr1 vs 1).
if tabix -l "$RAW" | grep -q '^chr'; then PFX="chr"; else PFX=""; fi
DEFAULT_CHROMS=""
for c in $(seq 1 22) X Y M; do DEFAULT_CHROMS+="${PFX}${c} "; done
CHROMS="${CHROMS:-$DEFAULT_CHROMS}"

present() { tabix -l "$RAW" | grep -qx "$1"; }

results=()
for chrom in $CHROMS; do
  if ! present "$chrom"; then
    echo "[skip] $chrom not in VCF"; continue
  fi
  echo "=== $chrom ==="
  raw_c="$OUT/${chrom}.raw.vcf.gz"
  ann_c="$OUT/${chrom}.ann.vcf.gz"
  bcftools view -r "$chrom" "$RAW" -Oz -o "$raw_c"
  tabix -f -p vcf "$raw_c"
  # annotate_vcf.sh does norm + chr<->Ensembl rename + SnpEff (SNPEFF_XMX) + vcfanno.
  "$HERE/annotate_vcf.sh" "$raw_c" "$REF" "$ann_c"
  # classify this chromosome on its own (peak RSS stays low: one chromosome's variants).
  python "$HERE/run_headless.py" "$ann_c" --hpo "$HPO" --out "$OUT" --sample-id "${chrom}"
  res="$OUT/${chrom}.ann_results.json"
  [[ -f "$res" ]] && results+=("$res") || echo "WARNING: no results.json for $chrom" >&2
  rm -f "$raw_c" "$raw_c".{tbi,csi}   # free disk as we go; keep the annotated per-chrom VCF
done

[[ ${#results[@]} -gt 0 ]] || { echo "ERROR: no per-chromosome results produced" >&2; exit 1; }

echo "=== merge ==="
python "$HERE/merge_wgs_results.py" "${results[@]}" --out "$OUT/merged.results.json"
# One whole-genome annotated VCF (chromosome order via bcftools concat, not glob order).
echo "=== concat annotated VCFs ==="
ls "$OUT"/*.ann.vcf.gz >/dev/null 2>&1 && \
  bcftools concat -Oz -o "$OUT/merged.annotated.vcf.gz" $(for c in $CHROMS; do echo "$OUT/${c}.ann.vcf.gz"; done 2>/dev/null | while read f; do [[ -f "$f" ]] && echo "$f"; done) && \
  tabix -f -p vcf "$OUT/merged.annotated.vcf.gz"

echo "Done -> $OUT/merged.results.json"
echo "Render the laudo from merged.results.json (same as any single-sample run)."
