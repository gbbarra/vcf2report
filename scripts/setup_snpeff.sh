#!/usr/bin/env bash
# Install SnpEff + the MANE GRCh38 database — the functional annotator for Stage 4.
#
# WHY this is needed at all: a raw caller's VCF (DRAGEN, GATK) carries only *quality*
# fields (AC/AF/DP/MQ/QD/FS/SOR). It has NO gene, consequence or HGVS — those come from
# a separate functional-annotation step. Without it the engine cannot evaluate
# PVS1/PM4/PP3 and the laudo has no c./p. notation. (The 1000G DRAGEN release publishes
# calls + QC only; it ships no annotation, so there is no shortcut.)
#
# WHY the MANE database specifically: the gnomAD store is sliced by MANE Select +
# MANE Plus Clinical (scripts/build_exome_bed.py, GENCODE v46). Annotating against the
# MANE db puts consequence and frequency on the SAME transcript — a transcript mismatch
# is a real source of PVS1 error. The .refseq flavour emits NM_ accessions, the
# convention clinical reports are written in.
#
#   bash scripts/setup_snpeff.sh
#
# Requires: java (brew install openjdk), curl. ~600 MB on disk, git-ignored.
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${VCF2REPORT_DATA:-$REPO/data}/tools"
JAR="$DEST/snpEff/snpEff.jar"
DB="${SNPEFF_DB:-GRCh38.mane.1.5.refseq}"

# Upstream moved hosts: the old snpeff.blob.core.windows.net is dead (NXDOMAIN).
# Do NOT substitute SourceForge or brewsci/bio — both are frozen at 4.3t (2017), whose
# newest human database is Ensembl 86 and which has no MANE databases at all.
URL="https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip"

command -v java >/dev/null 2>&1 || { echo "ERROR: java not found (brew install openjdk)"; exit 1; }

if [[ -f "$JAR" ]]; then
  echo "==> SnpEff already present: $JAR"
else
  mkdir -p "$DEST"
  echo "==> Downloading SnpEff core (~64 MB) to $DEST"
  curl -fL --retry 3 -o "$DEST/snpEff_core.zip" "$URL"
  unzip -q -o "$DEST/snpEff_core.zip" -d "$DEST"
  rm -f "$DEST/snpEff_core.zip"
fi

echo "==> SnpEff version: $(java -jar "$JAR" -version 2>&1 | head -1)"

if [[ -d "$DEST/snpEff/data/$DB" ]]; then
  echo "==> Database already present: $DB"
else
  echo "==> Downloading database $DB (~370 MB)"
  (cd "$DEST/snpEff" && java -jar snpEff.jar download -v "$DB" >/dev/null 2>&1)
fi

[[ -d "$DEST/snpEff/data/$DB" ]] || { echo "ERROR: database $DB failed to install"; exit 1; }

echo "==> Done: $JAR (db: $DB)"
echo "    Annotate with: bash scripts/annotate_vcf.sh RAW.vcf.gz OUT.annotated.vcf.gz"
