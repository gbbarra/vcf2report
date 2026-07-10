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

# SnpEff can be a PATH command (conda) OR a downloaded jar. Set SNPEFF_JAR=/path/
# snpEff.jar to use the zip install without any wrapper on PATH.
if [[ -n "${SNPEFF_JAR:-}" ]]; then
  SNPEFF=(java -Xmx8g -jar "$SNPEFF_JAR")
else
  SNPEFF=(snpEff)
fi

# RENAME_CHR=1 strips the "chr" prefix before SnpEff (DRAGEN VCFs are chr1.., the
# Ensembl SnpEff DB is 1..; mismatched names => SnpEff annotates nothing).
_norm_out="$RAW"

for tool in bcftools bgzip tabix vcfanno; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done
if [[ -z "${SNPEFF_JAR:-}" ]]; then
  command -v snpEff >/dev/null 2>&1 || { echo "ERROR: 'snpEff' not on PATH (or set SNPEFF_JAR=/path/snpEff.jar)" >&2; exit 1; }
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

echo "[1/3] Normalizing (split multiallelics + left-align) ..."
IN="$RAW"
if [[ "${RENAME_CHR:-0}" == "1" ]]; then
  echo "      RENAME_CHR=1 -> stripping 'chr' prefix to match the Ensembl SnpEff DB"
  bcftools annotate --rename-chrs \
    <(for c in $(seq 1 22) X Y MT; do echo "chr$c $c"; done; echo "chrM MT") \
    "$RAW" -Oz -o "$tmp/renamed.vcf.gz"
  IN="$tmp/renamed.vcf.gz"
fi
bcftools norm -m -any -f "$REF" --threads "$THREADS" -Oz -o "$tmp/norm.vcf.gz" "$IN"
bcftools index -t "$tmp/norm.vcf.gz"

echo "[2/3] SnpEff consequence + HGVS ($SNPEFF_DB) ..."
# -hgvs (default) emits HGVS.c/HGVS.p; -canon restricts to canonical transcripts.
"${SNPEFF[@]}" -noStats -hgvs "$SNPEFF_DB" "$tmp/norm.vcf.gz" | bgzip > "$tmp/snpeff.vcf.gz"
bcftools index -t "$tmp/snpeff.vcf.gz"

echo "[3/3] vcfanno population + clinical annotations ..."
# Write to a temp then move into place, so a mid-pipeline failure never leaves a
# truncated $OUT or clobbers a previously-good one.
vcfanno -p "$THREADS" "$VCFANNO_CONF" "$tmp/snpeff.vcf.gz" | bgzip > "$tmp/out.vcf.gz"
bcftools index -t "$tmp/out.vcf.gz"
mv -f "$tmp/out.vcf.gz" "$OUT"
mv -f "$tmp/out.vcf.gz.tbi" "$OUT.tbi"

echo "Done -> $OUT"
echo "Run: python scripts/run_headless.py $OUT --hpo patient_hpo_terms.txt"
