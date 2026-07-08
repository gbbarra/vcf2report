#!/usr/bin/env bash
# Annotate a raw exome VCF for vcf2report using a fully open-source, MIT/permissive,
# LOCAL toolchain (no per-variant network calls):
#
#   bcftools norm   (MIT)  -> normalize: split multiallelics + left-align indels
#   SnpEff          (MIT)  -> molecular consequence + HGVS (adds INFO/ANN)
#   vcfanno         (MIT)  -> gnomAD AF, ClinVar, REVEL/CADD from local files
#
# The result is an annotated VCF whose INFO vcf2report reads directly (see
# config.INFO_ALIASES), so the pipeline runs offline and fast on a real exome.
#
# Usage:
#   scripts/annotate_vcf.sh RAW.vcf.gz REF_GRCh38.fa OUT.annotated.vcf.gz
#
# Requirements (install once; all local): bcftools, snpEff, vcfanno, and the
# data files referenced by vcfanno.conf.toml (gnomAD sites VCF, ClinVar VCF, ...).
set -euo pipefail

RAW="${1:?raw VCF (bgzipped) required}"
REF="${2:?GRCh38 reference FASTA required}"
OUT="${3:?output path required}"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.105}"
VCFANNO_CONF="${VCFANNO_CONF:-$(dirname "$0")/vcfanno.conf.toml}"
THREADS="${THREADS:-4}"

for tool in bcftools vcfanno snpEff; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "[1/3] Normalizing (split multiallelics + left-align) ..."
bcftools norm -m -any -f "$REF" --threads "$THREADS" -Oz -o "$tmp/norm.vcf.gz" "$RAW"
bcftools index -t "$tmp/norm.vcf.gz"

echo "[2/3] SnpEff consequence + HGVS ($SNPEFF_DB) ..."
# -hgvs (default) emits HGVS.c/HGVS.p; -canon restricts to canonical transcripts.
snpEff -noStats -hgvs "$SNPEFF_DB" "$tmp/norm.vcf.gz" | bgzip > "$tmp/snpeff.vcf.gz"
bcftools index -t "$tmp/snpeff.vcf.gz"

echo "[3/3] vcfanno population + clinical annotations ..."
vcfanno -p "$THREADS" "$VCFANNO_CONF" "$tmp/snpeff.vcf.gz" | bgzip > "$OUT"
bcftools index -t "$OUT"

echo "Done -> $OUT"
echo "Run: python scripts/run_headless.py $OUT --hpo patient_hpo_terms.txt"
