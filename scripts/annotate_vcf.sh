#!/usr/bin/env bash
# Add functional annotation (gene + consequence + HGVS) to a raw exome VCF — Stage 4.
#
# A raw caller's VCF carries only quality fields (AC/AF/DP/MQ/QD/...). This adds the
# biology the engine needs: INFO/ANN with gene, molecular consequence and HGVS c./p.
# Without it, PVS1/PM4/PP3 are not evaluable and the laudo has no c./p. notation.
#
#   bcftools norm  (MIT) -> split multiallelics (+ left-align indels when a REF is given)
#   SnpEff         (MIT) -> consequence + HGVS on MANE transcripts (adds INFO/ANN)
#
# Population/clinical frequencies are NOT added here: the engine reads gnomAD,
# AlphaMissense and ClinVar straight from the Parquet stores (docs/DATA_ARCHITECTURE.md),
# so a vcfanno pass would be redundant. It still runs if vcfanno.conf.toml's data files
# happen to be present locally (VCFANNO=1 to force, =0 to skip).
#
# Usage:
#   scripts/annotate_vcf.sh RAW.vcf.gz OUT.annotated.vcf.gz [REF_GRCh38.fa]
#
# REF is optional: it only enables indel left-alignment in bcftools norm. Skipping it
# is safe for an already-normalized callset (the 1000G DRAGEN VCFs are).
#
# Install the annotator once with: bash scripts/setup_snpeff.sh
set -euo pipefail

RAW="${1:?raw VCF (bgzipped) required}"
OUT="${2:?output path required}"
REF="${3:-}"
REPO="$(cd "$(dirname "$0")/.." && pwd)"
SNPEFF_DB="${SNPEFF_DB:-GRCh38.mane.1.5.refseq}"
VCFANNO_CONF="${VCFANNO_CONF:-$REPO/scripts/vcfanno.conf.toml}"
THREADS="${THREADS:-4}"

# SnpEff resolution order: explicit SNPEFF_JAR, the setup_snpeff.sh location, then PATH.
DEFAULT_JAR="${VCF2REPORT_DATA:-$REPO/data}/tools/snpEff/snpEff.jar"
if [[ -n "${SNPEFF_JAR:-}" ]]; then
  SNPEFF=(java -Xmx8g -jar "$SNPEFF_JAR")
elif [[ -f "$DEFAULT_JAR" ]]; then
  SNPEFF=(java -Xmx8g -jar "$DEFAULT_JAR")
elif command -v snpEff >/dev/null 2>&1; then
  SNPEFF=(snpEff)
else
  echo "ERROR: SnpEff not found. Run: bash scripts/setup_snpeff.sh" >&2
  exit 1
fi

for tool in bcftools bgzip tabix; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' not found on PATH" >&2; exit 1; }
done

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# Chromosome naming must match or SnpEff silently annotates NOTHING (every record becomes
# ERROR_CHROMOSOME_NOT_FOUND). The GRCh38 MANE databases are Ensembl-style ("1"); DRAGEN and
# most clinical VCFs are UCSC-style ("chr1"). Detect both sides and rename in AND back out,
# so $OUT keeps the caller's original naming and still matches the chr-prefixed stores.
# set +o pipefail is scoped to this subshell: `head -1` closes the pipe, bcftools takes
# SIGPIPE (141), and pipefail would otherwise turn reading one line into a fatal error.
vcf_chr="$(set +o pipefail; bcftools view -H "$RAW" 2>/dev/null | head -1 | cut -f1)"
db_chr="$("${SNPEFF[@]}" genes2bed "$SNPEFF_DB" TP53 2>/dev/null | sed -n '2p' | cut -f1)"
vcf_ucsc=0; [[ "$vcf_chr" == chr* ]] && vcf_ucsc=1
db_ucsc=0;  [[ "$db_chr"  == chr* ]] && db_ucsc=1
if [[ -z "$db_chr" ]]; then
  echo "WARNING: could not probe '$SNPEFF_DB' chromosome naming (TP53 lookup failed);" >&2
  echo "         assuming Ensembl-style. Override with SNPEFF_DB_UCSC=1 if wrong." >&2
  db_ucsc="${SNPEFF_DB_UCSC:-0}"
fi

rename_needed=0
[[ "$vcf_ucsc" != "$db_ucsc" ]] && rename_needed=1

if [[ "$rename_needed" == 1 ]]; then
  for c in $(seq 1 22) X Y; do echo "chr$c $c"; done > "$tmp/to_ens.txt"; echo "chrM MT" >> "$tmp/to_ens.txt"
  for c in $(seq 1 22) X Y; do echo "$c chr$c"; done > "$tmp/to_ucsc.txt"; echo "MT chrM" >> "$tmp/to_ucsc.txt"
  if [[ "$vcf_ucsc" == 1 ]]; then IN_MAP="$tmp/to_ens.txt"; OUT_MAP="$tmp/to_ucsc.txt"
  else IN_MAP="$tmp/to_ucsc.txt"; OUT_MAP="$tmp/to_ens.txt"; fi
fi

echo "[1/3] Normalizing (split multiallelics${REF:+ + left-align}) ..."
norm_args=(-m -any --threads "$THREADS")
[[ -n "$REF" ]] && norm_args+=(-f "$REF")
bcftools norm "${norm_args[@]}" -Oz -o "$tmp/norm.vcf.gz" "$RAW"

echo "[2/3] SnpEff consequence + HGVS ($SNPEFF_DB) ..."
if [[ "$rename_needed" == 1 ]]; then
  echo "      chromosome naming: VCF=$vcf_chr db=$db_chr -> renaming in and back out"
  bcftools annotate --rename-chrs "$IN_MAP" "$tmp/norm.vcf.gz" -Oz -o "$tmp/in.vcf.gz"
else
  mv "$tmp/norm.vcf.gz" "$tmp/in.vcf.gz"
fi
# No -canon: with a MANE database it would keep only MANE Select and drop MANE Plus
# Clinical — the extra transcript that exists precisely for genes whose clinically
# relevant variants are not on Select. Dropping it would be a false negative. SnpEff
# orders ANN by severity, so the first entry stays the most severe consequence.
"${SNPEFF[@]}" -noStats -hgvs "$SNPEFF_DB" "$tmp/in.vcf.gz" > "$tmp/snpeff.vcf"
if [[ "$rename_needed" == 1 ]]; then
  bcftools annotate --rename-chrs "$OUT_MAP" "$tmp/snpeff.vcf" -Oz -o "$tmp/ann.vcf.gz"
else
  bgzip -c "$tmp/snpeff.vcf" > "$tmp/ann.vcf.gz"
fi

# Guard against the silent-failure mode: a naming or database mismatch yields a VCF that
# looks fine but carries no usable annotation. Fail loudly instead of shipping a laudo
# with no HGVS. (Records with ALT=* are spanning-deletion placeholders, not real calls,
# and are legitimately unannotatable — so require a bulk majority, not 100%.)
total="$(bcftools view -H "$tmp/ann.vcf.gz" | wc -l | tr -d ' ')"
with_ann="$(bcftools view -H "$tmp/ann.vcf.gz" | grep -c 'ANN=' || true)"
if [[ "$total" -gt 0 ]] && [[ "$with_ann" -lt $((total / 2)) ]]; then
  echo "ERROR: only $with_ann/$total records got ANN — annotation effectively failed." >&2
  echo "       Usual cause: chromosome naming or database mismatch (VCF=$vcf_chr db=$db_chr)." >&2
  exit 1
fi
echo "      annotated $with_ann/$total records"

echo "[3/3] vcfanno population + clinical annotations ..."
run_vcfanno="${VCFANNO:-auto}"
if [[ "$run_vcfanno" == "auto" ]]; then
  run_vcfanno=0
  if command -v vcfanno >/dev/null 2>&1 && [[ -f "$VCFANNO_CONF" ]]; then
    run_vcfanno=1
    while read -r f; do [[ -f "$f" ]] || run_vcfanno=0; done \
      < <(grep -E '^\s*file\s*=' "$VCFANNO_CONF" | sed -E 's/.*"(.*)".*/\1/')
  fi
fi
if [[ "$run_vcfanno" == 1 ]]; then
  vcfanno -p "$THREADS" "$VCFANNO_CONF" "$tmp/ann.vcf.gz" | bgzip > "$tmp/out.vcf.gz"
else
  echo "      skipped — the engine reads gnomAD/AlphaMissense/ClinVar from the Parquet stores"
  mv "$tmp/ann.vcf.gz" "$tmp/out.vcf.gz"
fi

# Write to a temp then move into place, so a mid-pipeline failure never leaves a
# truncated $OUT or clobbers a previously-good one.
bcftools index -t "$tmp/out.vcf.gz"
mkdir -p "$(dirname "$OUT")"
mv -f "$tmp/out.vcf.gz" "$OUT"
mv -f "$tmp/out.vcf.gz.tbi" "$OUT.tbi"

echo "Done -> $OUT"
echo "Run: python3 scripts/run_headless.py $OUT --hpo patient_hpo_terms.txt"
