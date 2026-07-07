"""Central configuration: paths, dataset locations, API endpoints, thresholds.

Everything that a deployment might want to point elsewhere lives here so the
rest of the code never hardcodes a path or a URL.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent.parent  # src/vcf2report -> repo root
DATA_DIR = Path(os.environ.get("VCF2REPORT_DATA", REPO_ROOT / "data"))
TEMPLATES_DIR = Path(os.environ.get("VCF2REPORT_TEMPLATES", REPO_ROOT / "templates"))
OUTPUT_DIR = Path(os.environ.get("VCF2REPORT_OUT", DATA_DIR / "out"))
CACHE_DIR = Path(os.environ.get("VCF2REPORT_CACHE", DATA_DIR / "cache"))

# Bundled local datasets (the offline fallback + secondary differentiators).
SAMPLE_VCF = DATA_DIR / "sample" / "sample_exome.vcf"
SAMPLE_HPO = DATA_DIR / "sample" / "patient_hpo_terms.txt"
CLINVAR_LOCAL = DATA_DIR / "clinvar" / "clinvar_grch38_slice.tsv"
GNOMAD_LOCAL = DATA_DIR / "gnomad" / "gnomad_cache.json"
ABRAOM_LOCAL = DATA_DIR / "abraom" / "abraom_sabe.tsv"
HPO_GENES_LOCAL = DATA_DIR / "hpo" / "genes_to_phenotype.tsv"
CONSTRAINT_LOCAL = DATA_DIR / "constraint" / "gene_constraint.tsv"
INSILICO_LOCAL = DATA_DIR / "insilico" / "insilico.tsv"

# ---------------------------------------------------------------------------
# Genome build — the whole pipeline assumes GRCh38 to match gnomAD r4 / ClinVar.
# ---------------------------------------------------------------------------
GENOME_BUILD = "GRCh38"

# ---------------------------------------------------------------------------
# Live API endpoints (used only when OFFLINE is false and deps are installed).
# ---------------------------------------------------------------------------
GNOMAD_API = "https://gnomad.broadinstitute.org/api"
GNOMAD_DATASET = "gnomad_r4"
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_API_KEY = os.environ.get("NCBI_API_KEY")  # optional; raises 3->10 req/s
NCBI_EMAIL = os.environ.get("NCBI_EMAIL", "vcf2report@example.org")
HPO_API = "https://ontology.jax.org/api"

# ---------------------------------------------------------------------------
# Behaviour flags
# ---------------------------------------------------------------------------
def offline() -> bool:
    """Cache-only mode. Set OFFLINE=1 for a network-independent demo."""
    return os.environ.get("OFFLINE", "").strip().lower() in {"1", "true", "yes"}


# ---------------------------------------------------------------------------
# Filtering / QC thresholds (documented so the report can cite them)
# ---------------------------------------------------------------------------
QC_MIN_DP = 10          # minimum read depth
QC_MIN_GQ = 20          # minimum genotype quality
QC_AB_MIN = 0.25        # het allele balance lower bound
QC_AB_MAX = 0.75        # het allele balance upper bound

# Population-frequency cutoffs for candidate rarity (popmax AF).
AF_DOMINANT_MAX = 0.001
AF_RECESSIVE_MAX = 0.005
# BA1 stand-alone benign threshold (Richards 2015).
AF_BA1 = 0.05
