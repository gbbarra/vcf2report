#!/usr/bin/env bash
# Publish the tell-free REALISTIC cohort as NEW assets on the existing SYN releases — added ALONGSIDE
# the current marked tarball (which is left untouched). Two assets per release: the raw tell-free
# VCFs and the same VCFs SnpEff-annotated. Truth stays external (planted_variants.tsv + per-sample
# sidecars). See data/synthetic_cohort/realistic/README.md.
#
#   scripts/publish_realistic_cohort.sh            # build + verify only (no upload)
#   scripts/publish_realistic_cohort.sh --upload   # also upload to GitHub
set -euo pipefail
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
R="$REPO_ROOT/data/synthetic_cohort/realistic"
UPLOAD=0; [ "${1:-}" = "--upload" ] && UPLOAD=1
command -v zstd >/dev/null || { echo "ERROR: zstd not found." >&2; exit 1; }
[ "$UPLOAD" = 1 ] && { command -v gh >/dev/null || { echo "ERROR: gh not found." >&2; exit 1; }; }
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

pack() {  # tag release
  local tag="$1" release="$2" prefix
  prefix="${release//-/_}"                       # syn-cohort-v2 -> syn_cohort_v2
  local dir="$R/$tag"
  local nraw ann
  nraw="$(find "$dir" -maxdepth 1 -name 'SYN-*.vcf.gz' | wc -l | tr -d ' ')"
  ann="$(find "$dir/annotated" -name 'SYN-*.annotated.vcf.gz' 2>/dev/null | wc -l | tr -d ' ')"
  [ "$nraw" -ge 100 ] && [ "$ann" -ge 100 ] || { echo "ERROR: $tag has raw=$nraw ann=$ann (<100)." >&2; exit 1; }

  echo ">>> $release: raw=$nraw annotated=$ann" >&2
  # raw tell-free + sidecars + manifest + README (exclude the annotated subdir)
  ( cd "$dir" && tar --exclude='._*' -cf - \
        $(ls SYN-*.vcf.gz SYN-*.vcf.gz.tbi SYN-*.planted.tsv SYN-*.hpo.txt) \
        -C "$R" README.md planted_variants.tsv ) \
    | zstd -3 -T0 -q -o "$WORK/${prefix}_realistic.tar.zst"
  # SnpEff-annotated realistic
  ( cd "$dir/annotated" && tar --exclude='._*' -cf - SYN-*.annotated.vcf.gz SYN-*.annotated.vcf.gz.tbi ) \
    | zstd -3 -T0 -q -o "$WORK/${prefix}_realistic_annotated.tar.zst"
  ( cd "$WORK" && shasum -a 256 "${prefix}_realistic.tar.zst" "${prefix}_realistic_annotated.tar.zst" > "${prefix}_realistic.SHA256" )
  echo "    $(du -h "$WORK/${prefix}_realistic.tar.zst" | cut -f1)  ${prefix}_realistic.tar.zst" >&2
  echo "    $(du -h "$WORK/${prefix}_realistic_annotated.tar.zst" | cut -f1)  ${prefix}_realistic_annotated.tar.zst" >&2
  # verify a sample tarball member is tell-free
  local probe
  probe="$(tar -tf <(zstd -dc "$WORK/${prefix}_realistic.tar.zst") | grep -m1 'SYN-.*\.vcf\.gz$')"
  echo "    sample member: $probe" >&2

  if [ "$UPLOAD" = 1 ]; then
    echo "    uploading to $release ..." >&2
    gh release upload "$release" \
      "$WORK/${prefix}_realistic.tar.zst" \
      "$WORK/${prefix}_realistic_annotated.tar.zst" \
      "$WORK/${prefix}_realistic.SHA256" --clobber
  fi
}

pack v1 syn-cohort-v2
pack v2 syn-cohort-expansion
echo "DONE ($([ "$UPLOAD" = 1 ] && echo uploaded || echo 'built + verified, NOT uploaded — pass --upload'))." >&2
