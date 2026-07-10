#!/usr/bin/env bash
# One-time setup: fetch the annotation databases vcf2report needs, LOCALLY.
# Run this ONCE per machine (not per analysis). After it, annotation is offline.
#
#   scripts/setup_data.sh [DATA_DIR] [--panel genes.bed]
#
# Downloads are large; use --panel to restrict gnomAD to a gene-panel BED (a few
# hundred MB instead of tens of GB) for a laptop / demo. Tools (bcftools, snpEff,
# vcfanno) are checked, not auto-installed — install via conda/bioconda:
#   conda install -c bioconda bcftools snpeff vcfanno htslib
set -euo pipefail

DATA_DIR="${1:-./annotation_data}"
PANEL_BED=""
[[ "${2:-}" == "--panel" ]] && PANEL_BED="${3:?panel BED path required after --panel}"
mkdir -p "$DATA_DIR"

echo "== Checking tools =="
missing=0
for t in bcftools tabix bgzip vcfanno snpEff; do
  if command -v "$t" >/dev/null 2>&1; then echo "  ok: $t"
  else echo "  MISSING: $t (conda install -c bioconda bcftools snpeff vcfanno htslib)"; missing=1; fi
done
[[ $missing -eq 1 ]] && echo "Install the missing tools, then re-run."

echo "== ClinVar (GRCh38, ~200 MB) =="
# -f/--fail so an HTTP error page isn't silently written as the .vcf.gz.
curl -fL -o "$DATA_DIR/clinvar_GRCh38.vcf.gz" \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
tabix -p vcf "$DATA_DIR/clinvar_GRCh38.vcf.gz" || true

echo "== SnpEff GRCh38 database =="
snpEff download -v GRCh38.105 || echo "  (run 'snpEff download GRCh38.105' manually if this failed)"

echo "== gnomAD v4 sites =="
if [[ -n "$PANEL_BED" ]]; then
  echo "  Panel mode: download gnomAD then subset to $PANEL_BED with:"
  echo "    bcftools view -R $PANEL_BED gnomad.sites.vcf.gz -Oz -o $DATA_DIR/gnomad_panel.vcf.gz"
  echo "  (Grab the sites VCF from https://gnomad.broadinstitute.org/downloads)"
else
  echo "  Full gnomAD sites VCF is tens of GB. Download from:"
  echo "    https://gnomad.broadinstitute.org/downloads  (v4 sites, GRCh38)"
  echo "  or re-run with:  scripts/setup_data.sh $DATA_DIR --panel your_genes.bed"
fi

echo "== ABraOM (Brazilian, SABE) =="
echo "  Download the release from http://abraom.ib.usp.br/ and convert to a"
echo "  bgzipped, tabix-indexed GRCh38 VCF at $DATA_DIR/abraom_GRCh38.vcf.gz"

echo ""
echo "Next: point scripts/vcfanno.conf.toml 'file' paths at $DATA_DIR, then:"
echo "  scripts/annotate_vcf.sh raw.vcf.gz GRCh38.fa out.annotated.vcf.gz"
