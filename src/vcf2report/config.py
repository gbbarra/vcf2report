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
HPO_GENES_LOCAL = DATA_DIR / "hpo" / "genes_to_phenotype.tsv.gz"
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
def _truthy(v: str | None) -> bool:
    return (v or "").strip().lower() in {"1", "true", "yes"}


def allow_network() -> bool:
    """Whether outbound calls to gnomAD/NCBI are permitted. OFF by default.

    Patient variant coordinates are sensitive genomic data, so egress is opt-IN:
    set ``VCF2REPORT_ALLOW_NETWORK=1`` to enable live lookups. ``OFFLINE=1`` always
    wins (forces no network) for a guaranteed-local demo.
    """
    if _truthy(os.environ.get("OFFLINE")):
        return False
    return _truthy(os.environ.get("VCF2REPORT_ALLOW_NETWORK"))


def offline() -> bool:
    """True when no network egress is allowed (the safe default)."""
    return not allow_network()


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

# ---------------------------------------------------------------------------
# ACMG SF v3.2 (Miller et al., 2023) secondary-findings genes. A P/LP variant in
# one of these, unrelated to the indication, is a reportable secondary finding
# (subject to the patient's opt-in). Curated subset of the ~81-gene list — verify
# and update against the current ACMG SF publication before clinical use.
# ---------------------------------------------------------------------------
ACMG_SF_GENES = {
    # Hereditary cancer
    "APC", "MUTYH", "BMPR1A", "SMAD4", "BRCA1", "BRCA2", "PALB2", "MLH1", "MSH2",
    "MSH6", "PMS2", "MEN1", "RET", "NF2", "SDHB", "SDHC", "SDHD", "SDHAF2", "MAX",
    "TMEM127", "VHL", "WT1", "TP53", "STK11", "PTEN", "CDH1", "RB1", "TSC1", "TSC2",
    # Cardiovascular
    "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYH11", "COL3A1", "LDLR", "APOB",
    "PCSK9", "MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "MYL3", "ACTC1", "PRKAG2",
    "MYL2", "LMNA", "RYR2", "PKP2", "DSP", "DSC2", "TMEM43", "DSG2", "KCNQ1", "KCNH2",
    "SCN5A", "CASQ2", "TRDN", "CALM1", "CALM2", "CALM3", "TNNC1", "BAG3", "DES", "FLNC",
    "RBM20", "TTN",
    # Malignant hyperthermia, metabolic, other
    "RYR1", "CACNA1S", "OTC", "GAA", "GLA", "ATP7B", "BTD", "RPE65", "TTR",
    "HFE", "ACVRL1", "ENG", "HNF1A",
}

# ---------------------------------------------------------------------------
# INFO field aliases for reading a *pre-annotated* VCF (SnpEff/VEP + vcfanno/
# bcftools). First matching key wins. Lets vcf2report consume a real annotated
# VCF fully offline — no per-variant DB lookups. Extend for your annotation.
# ---------------------------------------------------------------------------
INFO_ALIASES = {
    "gnomad_af": ["gnomad_AF", "gnomAD_AF", "gnomADg_AF", "gnomad4_AF", "AF_gnomad", "gnomad_af"],
    "gnomad_ac": ["gnomad_AC", "gnomAD_AC", "gnomad_ac"],
    "gnomad_an": ["gnomad_AN", "gnomAD_AN", "gnomad_an"],
    "gnomad_hom": ["gnomad_nhomalt", "gnomAD_nhomalt", "gnomad_hom", "nhomalt"],
    "abraom_af": ["ABraOM_AF", "abraom_AF", "ABRAOM_AF", "abraom_af"],
    "clinvar_sig": ["CLNSIG", "clinvar_CLNSIG", "clinvar_sig"],
    "clinvar_review": ["CLNREVSTAT", "clinvar_CLNREVSTAT"],
    "clinvar_disease": ["CLNDN", "clinvar_CLNDN"],
    "clinvar_accession": ["CLNVI", "ALLELEID", "clinvar_VCV"],
    "revel": ["REVEL", "dbNSFP_REVEL_score", "revel"],
    "cadd": ["CADD_PHRED", "CADD_phred", "cadd_phred", "CADD_PHRED_score"],
}
