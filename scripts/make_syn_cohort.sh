#!/usr/bin/env bash
# Build the SYN-001..N validation cohort: N DISTINCT 1000G DRAGEN v4.4.7 samples, each subset to the
# vendor-neutral MANE/GENCODE exome BED (data/gnomad/exome_hg38.bed — the same region the engine
# covers) and spiked with ONE distinct pathogenic ClinVar variant whose gene comes from a real GA4GH
# phenopacket case (carrying that case's HPO). No repeated background, no repeated variant — diverse
# real cases for validation. Writes SYN-00N.synthetic.vcf.gz + SYN-00N.hpo.txt + truth.tsv.
#
# SELF-CONTAINED: needs only curl, bcftools, bgzip, tabix, python3 (no aws CLI, no reference FASTA).
# The DRAGEN VCFs are streamed from the public S3 bucket over HTTPS; the ClinVar VCF is auto-downloaded.
#
#   bash scripts/make_syn_cohort.sh                 # all 100
#   N=3 bash scripts/make_syn_cohort.sh             # just the first 3 (a quick smoke test)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
COHORT="${1:-$REPO/data/synthetic_cohort/cohort.tsv}"
OUT="${2:-$REPO/data/synthetic_cohort}"
BED="$REPO/data/gnomad/exome_hg38.bed"
BUCKET="1000genomes-dragen-v4-4-7"
PIPE="data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"
LIMIT="${N:-0}"   # 0 = all

for t in curl bcftools bgzip tabix python3; do
  command -v "$t" >/dev/null || { echo "ERROR: '$t' not on PATH" >&2; exit 1; }
done
[ -f "$BED" ]    || { echo "ERROR: MANE BED not found: $BED" >&2; exit 1; }
[ -f "$COHORT" ] || { echo "ERROR: cohort not found: $COHORT" >&2; exit 1; }
mkdir -p "$OUT"

# ClinVar GRCh38 VCF (for the spike's CLNSIG/CLNREVSTAT/CLNDN). Use $CLINVAR_VCF, else the repo copy,
# else download it once (~180 MB, public domain — NCBI).
CLINVAR="${CLINVAR_VCF:-$REPO/data/clinvar/clinvar.vcf.gz}"
if [ ! -f "$CLINVAR" ]; then
  echo "== fetching ClinVar GRCh38 VCF (~180 MB, one-time) -> $CLINVAR"
  mkdir -p "$(dirname "$CLINVAR")"
  curl -fL --retry 3 -o "$CLINVAR" "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
fi

TRUTH="$OUT/truth.tsv"; [ -f "$TRUTH" ] || echo -e "syn_id\tsample\tgene\tvariant\tconsequence\thpo_n\texome_variants\tdisease" > "$TRUTH"

count=0
tail -n +2 "$COHORT" | while IFS=$'\t' read -r syn sample gene chrom pos ref alt cons disease hpo; do
  [ -n "${syn:-}" ] || continue
  count=$((count+1)); [ "$LIMIT" -gt 0 ] && [ "$count" -gt "$LIMIT" ] && break
  if [ -f "$OUT/$syn.synthetic.vcf.gz" ]; then echo "==== [$syn] present — skip"; continue; fi
  echo "==== [$syn] $sample  $gene  $chrom:$pos $ref>$alt ($cons) ===="
  work="$OUT/.work/$syn"; mkdir -p "$work"

  # 1) DRAGEN VCF over HTTPS — deterministic key, else list to find it
  url="https://$BUCKET.s3.amazonaws.com/$PIPE/$sample/$sample.hard-filtered.vcf.gz"
  if ! curl -fsI --max-time 30 "$url" >/dev/null 2>&1; then
    key="$(curl -s --max-time 30 "https://$BUCKET.s3.amazonaws.com/?list-type=2&prefix=$PIPE/$sample/&max-keys=200" \
           | grep -oE "<Key>[^<]+hard-filtered\.vcf\.gz</Key>" | sed 's/<[^>]*>//g' | head -1 || true)"
    [ -n "$key" ] || { echo "  WARN: no DRAGEN VCF for $sample — skipping"; rm -rf "$work"; continue; }
    url="https://$BUCKET.s3.amazonaws.com/$key"
  fi
  curl -fL --retry 3 --max-time 1800 -o "$work/raw.vcf.gz" "$url"

  # 2) split multiallelics (no FASTA needed) + subset to the MANE BED — one streaming pass
  bcftools norm -m -any "$work/raw.vcf.gz" -Ou 2>/dev/null \
    | bcftools view -T "$BED" -Oz -o "$work/exome.vcf.gz"

  # 3) spike the gene's pathogenic ClinVar variant (carries CLNSIG/CLNREVSTAT/CLNDN -> PP5 + disease)
  printf 'gene\tcategory\tzygosity\n%s\tprimary\thet\n' "$gene" > "$work/targets.tsv"
  python3 "$REPO/scripts/spike_pathogenic.py" \
    --exome "$work/exome.vcf.gz" --clinvar "$CLINVAR" --targets "$work/targets.tsv" \
    --sample-id "$syn" --out "$OUT/$syn.synthetic.vcf"
  bgzip -f "$OUT/$syn.synthetic.vcf"; tabix -f -p vcf "$OUT/$syn.synthetic.vcf.gz"

  # 4) HPO + truth
  printf '%s\n' "$hpo" | tr ',' '\n' | sed '/^$/d' > "$OUT/$syn.hpo.txt"
  vc="$(bgzip -dc "$OUT/$syn.synthetic.vcf.gz" | grep -c '^[^#]' || echo 0)"
  hn="$(wc -l < "$OUT/$syn.hpo.txt" | tr -d ' ')"
  echo -e "$syn\t$sample\t$gene\t$chrom:$pos:$ref:$alt\t$cons\t$hn\t$vc\t$disease" >> "$TRUTH"
  echo "  -> $OUT/$syn.synthetic.vcf.gz  ($vc variants, $hn HPO)"
  rm -rf "$work"
done
rmdir "$OUT/.work" 2>/dev/null || true

echo ""
echo "Done -> $OUT  (ground truth: $TRUTH)"
echo "Validate: bash scripts/validate_cohort.sh"
