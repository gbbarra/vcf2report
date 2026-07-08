#!/usr/bin/env bash
# Build N synthetic case exomes from public 1000G DRAGEN samples:
#   download per-sample WGS VCF -> normalize -> subset to IDT xGen Exome v2 ->
#   spike real ClinVar pathogenic variants -> de-identify -> bgzip+tabix.
#
# Real, diverse background (5 super-populations) + planted pathogenic variants,
# so the demo reliably shows Pathogenic/Likely-Pathogenic calls on realistic data.
#
# Requirements (all local, run on YOUR machine — this needs AWS S3 + the IDT BED,
# which the sandbox can't reach): awscli, bcftools, bgzip, tabix, python3.
#
#   ./make_synthetic_exomes.sh
#
# The per-sample spike targets + HPO are read from data/synthetic/SYN-00N.* (the
# validated definitions: primary DEE gene is LoF-intolerant + phenotype-matching,
# secondary ACMG SF gene is LoF-intolerant + has zero seizure-HPO overlap so it
# routes to "secondary" and not "primary" under the full real HPO). See
# data/synthetic/README.md.
set -euo pipefail

# ---------------------------------------------------------------------------
# CONFIG — edit these.
# ---------------------------------------------------------------------------
# 1000G DRAGEN 4.4.7 reanalysis — public AWS Open Data bucket (no credentials).
# The per-sample small-variant VCF key layout below is VERIFIED (NA12878 etc. return
# HTTP 200). The bucket name has no dots, so the HTTPS virtual-hosted form also works
# with plain curl if you don't have the aws CLI. Confirm a sample exists with:
#   aws s3 ls --no-sign-request s3://1000genomes-dragen-v4-4-7/data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2/NA12878/
S3_BUCKET="s3://1000genomes-dragen-v4-4-7"
S3_PIPE="data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2"
S3_VCF_TEMPLATE="${S3_BUCKET}/${S3_PIPE}/{SAMPLE}/{SAMPLE}.hard-filtered.vcf.gz"
# (older fallback, GRCh38, per-sample: s3://1000genomes-dragen/data/dragen-3.5.7b/hg38_altaware_nohla-cnv-anchored/{SAMPLE}/{SAMPLE}.hard-filtered.vcf.gz
#  — but NA19240 is ABSENT there; use NA19238/NA19239 (YRI) if you downgrade.)

# IDT xGen Exome Hyb Panel v2 targets (hg38). A byte-exact public mirror of IDT's
# official BED (verified reachable); or download from IDT and set a local path.
IDT_BED="idt_xgen_exome_v2_hg38.bed"
IDT_BED_URL="https://raw.githubusercontent.com/icgc-argo-workflows/target_seq_bed_files/main/TEST-INTL/WXS/xGenTM_Exome_Hyb_panel_v2/primary_target_regions.bed"
REF_FASTA="GRCh38.fa"                     # for bcftools norm (indexed .fai)
# Reachable GRCh38 FASTA mirror: https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta (+ .fai)
CLINVAR_VCF="clinvar_GRCh38.vcf.gz"       # ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
# Reachable ClinVar GRCh38 mirror (uncompressed .vcf; bgzip+tabix it locally for vcfanno):
#   https://storage.googleapis.com/gcp-public-data--broad-references/funcotator/funcotator_dataSources.v1.8.hg38.20230908g/clinvar/hg38/clinvar_20230717_hg38.vcf
HERE="$(cd "$(dirname "$0")" && pwd)"
SPIKE_SCRIPT="$HERE/spike_pathogenic.py"
SYNDIR="$HERE/../data/synthetic"          # validated targets + HPO live here
OUTDIR="synthetic_exomes"

# sample : SYN id : super-population label. The primary/secondary genes come from
# data/synthetic/${SYNID}.targets.tsv (validated), not from this array.
# NOTE: NA12878 (CEU) is confirmed 1000G + GIAB. VERIFY the others exist in the
# bucket (aws s3 ls ...) and swap freely; the population labels are the intent.
SAMPLES=(
  "NA12878:SYN-001:EUR (CEU) — has GIAB truth set"
  "NA19240:SYN-002:AFR (YRI)"
  "NA18939:SYN-003:EAS (JPT)"
  "NA20845:SYN-004:SAS (GIH)"
  "HG01112:SYN-005:AMR (CLM) — closest to admixed Brazilian ancestry"
)
# ---------------------------------------------------------------------------

mkdir -p "$OUTDIR"
for tool in aws bcftools bgzip tabix python3; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not on PATH" >&2; exit 1; }
done
if [[ ! -f "$IDT_BED" ]]; then
  echo "IDT BED not found locally; fetching mirror -> $IDT_BED"
  curl -fsSL "$IDT_BED_URL" -o "$IDT_BED" || { echo "ERROR: could not fetch IDT BED" >&2; exit 1; }
fi
[[ -f "$CLINVAR_VCF" ]] || { echo "ERROR: ClinVar VCF not found: $CLINVAR_VCF" >&2; exit 1; }

i=0
for entry in "${SAMPLES[@]}"; do
  i=$((i+1))
  SAMPLE="${entry%%:*}"; rest="${entry#*:}"
  SYNID="${rest%%:*}"; POP="${rest#*:}"
  TARGETS="$SYNDIR/${SYNID}.targets.tsv"
  HPO_SRC="$SYNDIR/${SYNID}.hpo.txt"
  [[ -f "$TARGETS" ]] || { echo "ERROR: missing $TARGETS" >&2; exit 1; }
  echo ""; echo "==================================================================="
  echo "[$i/${#SAMPLES[@]}] $SAMPLE  ->  $SYNID   [$POP]"
  sed 's/^/      target: /' "$TARGETS"
  echo "==================================================================="
  work="$OUTDIR/$SYNID"; mkdir -p "$work"

  s3uri="${S3_VCF_TEMPLATE//\{SAMPLE\}/$SAMPLE}"
  echo "[1/4] download  $s3uri"
  aws s3 cp --no-sign-request "$s3uri" "$work/raw.vcf.gz"
  tabix -f -p vcf "$work/raw.vcf.gz" 2>/dev/null || bcftools index -t "$work/raw.vcf.gz"

  echo "[2/4] normalize + subset to IDT exome v2"
  # -T (targets, streaming) not -R (regions, needs an index) because the input is
  # a pipe. bcftools accepts a BED for -T; overlap-filters as it streams.
  bcftools norm -m -any -f "$REF_FASTA" -Ou "$work/raw.vcf.gz" \
    | bcftools view -T "$IDT_BED" -Oz -o "$work/exome.vcf.gz"
  tabix -f -p vcf "$work/exome.vcf.gz"
  echo "      exome variants: $(bcftools view -H "$work/exome.vcf.gz" | wc -l)"

  echo "[3/4] spike pathogenic (targets from $TARGETS)"
  python3 "$SPIKE_SCRIPT" \
    --exome "$work/exome.vcf.gz" --clinvar "$CLINVAR_VCF" \
    --targets "$TARGETS" --sample-id "$SYNID" \
    --out "$work/${SYNID}.synthetic.vcf"

  echo "[4/4] sort + bgzip + index"
  bcftools sort "$work/${SYNID}.synthetic.vcf" -Oz -o "$OUTDIR/${SYNID}.synthetic.vcf.gz"
  tabix -f -p vcf "$OUTDIR/${SYNID}.synthetic.vcf.gz"
  cp "$HPO_SRC" "$OUTDIR/${SYNID}.hpo.txt"    # validated seizure-specific HPO
  rm -rf "$work"
  echo "  -> $OUTDIR/${SYNID}.synthetic.vcf.gz  (+ ${SYNID}.hpo.txt)"
done

echo ""; echo "Done. ${#SAMPLES[@]} synthetic exomes in $OUTDIR/"
echo "NEXT: annotate each (SnpEff/vcfanno) then run vcf2report, e.g.:"
echo "  scripts/annotate_vcf.sh $OUTDIR/SYN-001.synthetic.vcf.gz $REF_FASTA SYN-001.annotated.vcf.gz"
echo "  python scripts/run_headless.py SYN-001.annotated.vcf.gz --hpo $OUTDIR/SYN-001.hpo.txt"
