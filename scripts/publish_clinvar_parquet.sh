#!/usr/bin/env bash
# Publish the locally-built ClinVar Parquet store as a checksummed GitHub Release asset,
# so anyone can fetch + verify it without downloading the ClinVar VCF and rebuilding. ClinVar
# releases WEEKLY — re-run this weekly (build_clinvar_local.py -> build_clinvar_parquet.py ->
# this) to refresh the asset IN PLACE (--clobber), keeping the fetch URL stable.
#
#   scripts/publish_clinvar_parquet.sh [STORE_DIR]   # default: data/clinvar/clinvar_parquet
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STORE="${1:-$REPO_ROOT/data/clinvar/clinvar_parquet}"
TAG="clinvar-parquet-latest"
ASSET="clinvar_parquet.tar.zst"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

[ -d "$STORE" ] || { echo "ERROR: store not found: $STORE — build it first (build_clinvar_local.py -> build_clinvar_parquet.py)" >&2; exit 1; }
[ -f "$STORE/_manifest.json" ] || { echo "ERROR: $STORE/_manifest.json missing — the build looks incomplete." >&2; exit 1; }
command -v gh   >/dev/null || { echo "ERROR: gh (GitHub CLI) not found — https://cli.github.com" >&2; exit 1; }
command -v zstd >/dev/null || { echo "ERROR: zstd not found (brew install zstd / apt install zstd)." >&2; exit 1; }

echo "Packaging $STORE -> $ASSET ..." >&2
# Parquet is already columnar-compressed; -3 is fast and the archive stays ~store size.
tar --exclude='._*' -C "$(dirname "$STORE")" -cf - "$(basename "$STORE")" | zstd -3 -T0 -o "$WORK/$ASSET" -q
( cd "$WORK" && shasum -a 256 "$ASSET" > SHA256SUMS )
echo "  size: $(du -h "$WORK/$ASSET" | cut -f1)" >&2
cat "$WORK/SHA256SUMS" >&2

BUILT="$(python3 -c "import json;print(json.load(open('$STORE/_manifest.json'))['built_at'])" 2>/dev/null || echo unknown)"
echo "Publishing release $TAG (built_at $BUILT) ..." >&2
if gh release view "$TAG" >/dev/null 2>&1; then
  gh release upload "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" --clobber
else
  gh release create "$TAG" "$WORK/$ASSET" "$WORK/SHA256SUMS" \
    --title "ClinVar GRCh38 Parquet (vcf2report)" \
    --notes "ClinVar GRCh38 variant-classification Parquet for vcf2report (public domain — NCBI). Rebuilt WEEKLY from ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz. Fetch + verify: scripts/fetch_clinvar_parquet.sh"
fi
echo "Done. Anyone can now: scripts/fetch_clinvar_parquet.sh" >&2
