#!/usr/bin/env bash
# Build the SYN-001..N validation cohort: N DISTINCT 1000G DRAGEN v4.4.7 samples, each subset to
# the vendor-neutral MANE/GENCODE exome BED (data/gnomad/exome_hg38.bed — the SAME region the engine
# covers) and spiked with ONE distinct pathogenic ClinVar variant whose gene comes from a real
# GA4GH-phenopacket case (carrying that case's HPO). No repeated background, no repeated variant —
# diverse real cases for validation. Writes SYN-00N.synthetic.vcf.gz + SYN-00N.hpo.txt + truth.tsv.
#
#   bash scripts/make_syn_cohort.sh [cohort.tsv] [OUT_DIR]
#
# Run on YOUR machine (the sandbox can't reach S3 fast enough). Prereqs: awscli, bcftools, bgzip,
# tabix, python3, a GRCh38 FASTA (REF_FASTA=…, indexed .fai), a ClinVar GRCh38 VCF (CLINVAR_VCF=…).
# cohort.tsv columns: syn_id  sample  gene  chrom  pos  ref  alt  consequence  disease  hpo(comma)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
COHORT="${1:-$REPO/data/synthetic_cohort/cohort.tsv}"
OUT="${2:-$REPO/data/synthetic_cohort}"
BED="$REPO/data/gnomad/exome_hg38.bed"
REF="${REF_FASTA:?set REF_FASTA=/path/GRCh38.fa (with an indexed .fai)}"
CLINVAR="${CLINVAR_VCF:?set CLINVAR_VCF=/path/clinvar_GRCh38.vcf.gz}"
BUCKET="1000genomes-dragen-v4-4-7"
PIPE="data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"

for t in aws bcftools bgzip tabix python3; do
  command -v "$t" >/dev/null || { echo "ERROR: '$t' not on PATH" >&2; exit 1; }
done
[ -f "$BED" ]     || { echo "ERROR: MANE BED not found: $BED" >&2; exit 1; }
[ -f "$COHORT" ]  || { echo "ERROR: cohort not found: $COHORT" >&2; exit 1; }
mkdir -p "$OUT"
TRUTH="$OUT/truth.tsv"; echo -e "syn_id\tsample\tgene\tvariant\tconsequence\thpo_n\texome_variants\tdisease" > "$TRUTH"

tail -n +2 "$COHORT" | while IFS=$'\t' read -r syn sample gene chrom pos ref alt cons disease hpo; do
  [ -n "${syn:-}" ] || continue
  if [ -f "$OUT/$syn.synthetic.vcf.gz" ]; then echo "==== [$syn] present — skip"; continue; fi
  echo "==== [$syn] $sample  $gene  $chrom:$pos $ref>$alt ($cons) ===="
  work="$OUT/.work/$syn"; mkdir -p "$work"

  # 1) locate + download the sample's DRAGEN hard-filtered VCF (+ its index)
  key="$(aws s3 ls --no-sign-request --recursive "s3://$BUCKET/$PIPE/$sample/" 2>/dev/null \
         | grep -oE "[^ ]+hard-filtered\.vcf\.gz$" | head -1 || true)"
  if [ -z "$key" ]; then echo "  WARN: no DRAGEN VCF for $sample — skipping"; rm -rf "$work"; continue; fi
  aws s3 cp --no-sign-request "s3://$BUCKET/$key"     "$work/raw.vcf.gz"     >/dev/null
  aws s3 cp --no-sign-request "s3://$BUCKET/$key.tbi" "$work/raw.vcf.gz.tbi" >/dev/null 2>&1 \
    || tabix -f -p vcf "$work/raw.vcf.gz"

  # 2) normalize + subset to the MANE/GENCODE exome BED (one streaming pass; vendor-neutral)
  bcftools norm -m -any -f "$REF" -Ou "$work/raw.vcf.gz" \
    | bcftools view -T "$BED" -Oz -o "$work/exome.vcf.gz"

  # 3) spike the gene's pathogenic ClinVar variant (carries CLNSIG/CLNREVSTAT/CLNDN -> PP5 + disease),
  #    de-identified to $syn
  printf 'gene\tcategory\tzygosity\n%s\tprimary\thet\n' "$gene" > "$work/targets.tsv"
  python3 "$REPO/scripts/spike_pathogenic.py" \
    --exome "$work/exome.vcf.gz" --clinvar "$CLINVAR" --targets "$work/targets.tsv" \
    --sample-id "$syn" --out "$OUT/$syn.synthetic.vcf"
  bgzip -f "$OUT/$syn.synthetic.vcf"; tabix -f -p vcf "$OUT/$syn.synthetic.vcf.gz"

  # 4) HPO file (the phenopacket case's terms) + truth row
  printf '%s\n' "$hpo" | tr ',' '\n' | sed '/^$/d' > "$OUT/$syn.hpo.txt"
  vc="$(bgzip -dc "$OUT/$syn.synthetic.vcf.gz" | grep -c '^[^#]' || echo 0)"
  hn="$(wc -l < "$OUT/$syn.hpo.txt" | tr -d ' ')"
  echo -e "$syn\t$sample\t$gene\t$chrom:$pos:$ref:$alt\t$cons\t$hn\t$vc\t$disease" >> "$TRUTH"
  echo "  -> $OUT/$syn.synthetic.vcf.gz  ($vc variants, $hn HPO)"
  rm -rf "$work"
done
rmdir "$OUT/.work" 2>/dev/null || true

echo ""
echo "Cohort built -> $OUT  (ground truth: $TRUTH)"
echo "Validate one:  python3 scripts/run_headless.py $OUT/SYN-001.synthetic.vcf.gz --hpo $OUT/SYN-001.hpo.txt"
echo "Validate all:  bash scripts/validate_cohort.sh   (surfaces the planted gene per case vs truth.tsv)"
