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
# Reduced local gnomAD frequency table (bgzipped TSV + .tbi), built by
# scripts/build_gnomad_local.py. Offline, authoritative when present. Overridable so
# a large full build can live on an external disk. Schema (tabix -s1 -b2 -e2):
#   #chrom  pos  ref  alt  af  ac  an  hom  faf95  pop
GNOMAD_LOCAL_TABIX = Path(os.environ.get(
    "VCF2REPORT_GNOMAD_TABIX", str(DATA_DIR / "gnomad" / "gnomad_freq.local.tsv.gz")))
# gnomAD frequencies as a DuckDB/Parquet store (built by scripts/build_gnomad_parquet.py,
# or an existing lakehouse gnomad_freq.parquet). A single .parquet file or a Hive-
# partitioned dir (chrom=chrN/). Whole-exome frequencies come from one vectorised join
# in ~seconds, offline. Absent -> feature off. Overridable env.
GNOMAD_PARQUET = os.environ.get("VCF2REPORT_GNOMAD_PARQUET") or None
ABRAOM_LOCAL = DATA_DIR / "abraom" / "abraom_sabe.tsv"
HPO_GENES_LOCAL = DATA_DIR / "hpo" / "genes_to_phenotype.tsv.gz"
HPO_GRAPH_LOCAL = DATA_DIR / "hpo" / "hpo_graph.tsv.gz"  # ontology + IC (build_hpo_graph.py)
# Report routing: a gene is "phenotype-related" (-> primary findings) when its
# SINGLE strongest patient<->gene term similarity (hpo_best_match) is at/above this.
# Routing on the best match, not the average, keeps a gene that strongly explains one
# key phenotype in primary instead of diluting it out on a phenotype-rich case. Set
# above the incidental ontology noise floor (unrelated genes peak ~0.4) and below a
# real moderate match; distinct from the PP4 evidence bar (0.6), which uses the average.
HPO_RELATED_MIN = 0.5
CONSTRAINT_LOCAL = DATA_DIR / "constraint" / "gene_constraint.tsv.gz"
INSILICO_LOCAL = DATA_DIR / "insilico" / "insilico.tsv"
# AlphaMissense hg38 predictions (CC BY 4.0) — tabix-indexed, fetched once via
# scripts/fetch_alphamissense.sh. Absent by default; the client degrades to None.
ALPHAMISSENSE_LOCAL = Path(os.environ.get(
    "VCF2REPORT_ALPHAMISSENSE", DATA_DIR / "alphamissense" / "AlphaMissense_hg38.tsv.gz"))

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
# ACMG combining model — Richards 2015 Table 5 (default) vs the ClinGen/Tavtigian
# naturally-scaled POINTS system (Tavtigian et al., Genet Med 2020) together with
# the ClinGen SVI 2020 refinement that downgrades PM2 from Moderate to Supporting.
# Toggle with VCF2REPORT_ACMG_MODEL=clingen. Deterministic either way.
# ---------------------------------------------------------------------------
def acmg_model() -> str:
    m = (os.environ.get("VCF2REPORT_ACMG_MODEL") or "richards").strip().lower()
    return "clingen" if m in ("clingen", "clingen2020", "points") else "richards"


def pm2_strength() -> str:
    """PM2 applied strength: Moderate (Richards) or Supporting (ClinGen SVI 2020)."""
    return "supporting" if acmg_model() == "clingen" else "moderate"


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
# BS1 — "allele frequency greater than expected for the disorder" (Richards 2015,
# refined by the filtering-AF framework of Whiffin et al., Genet Med 2017). The
# threshold is disorder-dependent, NOT a single universal number: a dominant
# condition tolerates far less allele frequency than a recessive one (where the
# carrier frequency is expected to be higher). We key the cutoff on the gene's
# mode of inheritance and fall back to a conservative default when it is unknown.
# These are pragmatic gnomAD-based heuristics, not a per-gene max-credible-AF
# derived from prevalence/penetrance — the report labels them as such.
# ---------------------------------------------------------------------------
BS1_AF_DOMINANT = 0.001     # rare dominant disorder: even ~0.1% is too common
BS1_AF_RECESSIVE = 0.01     # recessive: carrier frequency runs higher
BS1_AF_DEFAULT = 0.005      # inheritance unknown → conservative middle ground

# Curated mode-of-inheritance for the genes exercised by the demo/secondary-finding
# cases (and a few common ones). "AD"/"AR"/"XL"; genes absent here resolve to the
# conservative BS1_AF_DEFAULT. A curated subset — extend/verify against a source
# such as OMIM/Genomics England PanelApp before clinical use.
GENE_INHERITANCE = {
    # DEE / seizure primaries (haploinsufficiency-driven dominant)
    "SCN1A": "AD", "SCN2A": "AD", "KCNQ2": "AD", "STXBP1": "AD",
    "SLC2A1": "AD", "CACNA1A": "AD", "PAX6": "AD",
    # ACMG SF secondaries used in the synthetic cases
    "RB1": "AD", "APC": "AD", "STK11": "AD", "WT1": "AD", "FBN1": "AD",
    # a couple of well-known recessive genes for contrast/tests
    "CFTR": "AR", "HFE": "AR", "GJB2": "AR", "MUTYH": "AR",
}


def bs1_af_cutoff(gene: str | None) -> tuple[float, str | None]:
    """BS1 allele-frequency cutoff and the mode of inheritance it came from.

    Returns ``(cutoff, moi)`` where ``moi`` is "AD"/"AR"/"XL" when the gene is in
    the curated map, else ``None`` (→ conservative default cutoff).
    """
    moi = GENE_INHERITANCE.get((gene or "").upper()) if gene else None
    if moi == "AR":
        return BS1_AF_RECESSIVE, moi
    if moi in ("AD", "XL"):
        return BS1_AF_DOMINANT, moi
    return BS1_AF_DEFAULT, None


# ---------------------------------------------------------------------------
# PM2 — "absent/ultra-rare in population databases" (ClinGen SVI: supporting).
# Like BS1, the credible ceiling is disorder-dependent: a recessive disorder
# tolerates a higher AF (carriers are common in the general population) before a
# variant stops looking rare, whereas a dominant disorder needs near-absence.
# ---------------------------------------------------------------------------
PM2_AF_DOMINANT = 1e-4      # dominant: absent / ultra-rare
PM2_AF_RECESSIVE = 1e-3     # recessive: carrier frequency tolerated
PM2_AF_DEFAULT = 1e-4       # inheritance unknown → strict default


def pm2_af_ceiling(gene: str | None) -> tuple[float, str | None]:
    """PM2 rarity ceiling and the mode of inheritance it came from.

    Returns ``(ceiling, moi)``; genes absent from the curated map resolve to the
    strict default so PM2 is never granted more liberally than we can justify.
    """
    moi = GENE_INHERITANCE.get((gene or "").upper()) if gene else None
    if moi == "AR":
        return PM2_AF_RECESSIVE, moi
    if moi in ("AD", "XL"):
        return PM2_AF_DOMINANT, moi
    return PM2_AF_DEFAULT, None

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
    # gnomAD filtering AF (95% CI, grpmax) — the ClinGen-recommended field for BS1/BA1.
    "gnomad_faf95": ["gnomad_faf95", "fafmax_faf95_max", "faf95_max", "faf95_grpmax",
                     "gnomad_faf95_max", "AF_grpmax_faf95", "faf95"],
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
    "am_pathogenicity": ["am_pathogenicity", "AlphaMissense", "AlphaMissense_score",
                         "dbNSFP_AlphaMissense_score", "alphamissense"],
    "am_class": ["am_class", "AlphaMissense_class", "AlphaMissense_pred", "am_classification"],
}

# ---------------------------------------------------------------------------
# AlphaMissense (am_pathogenicity in 0..1) -> ACMG/AMP PP3/BP4 evidence strength.
#
# ClinGen's 2024 recalibration (Schmidt et al., Genet Med 2025) found AlphaMissense
# can reach the STRONG level for pathogenicity (PP3) and the MODERATE level for
# benignity (BP4) — but only at score cutoffs MORE stringent than the tool's own
# 0.564 / 0.34 class boundaries. The exact per-strength cutoffs from that paper are
# not reproduced here; the values below are documented SEED thresholds. They are
# meant to be calibrated empirically against the concordance panel (raise them
# until gross discordances stay at zero) and VERIFIED against the ClinGen table
# before any clinical use. Native AlphaMissense classes: <0.34 likely_benign,
# 0.34-0.564 ambiguous, >0.564 likely_pathogenic.
# ---------------------------------------------------------------------------
AM_PP3_STRONG = 0.99        # >= this -> PP3 at Strong (with PM2 -> Likely Pathogenic)
AM_PP3_MODERATE = 0.90      # >= this -> PP3 at Moderate
AM_PP3_SUPPORTING = 0.564   # >= this -> PP3 at Supporting (tool's likely_pathogenic)
# Richards 2015 Table 5 has no benign "moderate" bucket, so AlphaMissense benign
# evidence is capped at Supporting here (a documented limitation of the classic
# combining rules vs. ClinGen's points framework).
AM_BP4_SUPPORTING = 0.34    # <= this -> BP4 at Supporting (tool's likely_benign)


def am_pp3_strength(am: float | None) -> str | None:
    """PP3 evidence strength for an AlphaMissense score, or None if it doesn't apply."""
    if am is None:
        return None
    if am >= AM_PP3_STRONG:
        return "strong"
    if am >= AM_PP3_MODERATE:
        return "moderate"
    if am >= AM_PP3_SUPPORTING:
        return "supporting"
    return None


def am_bp4_strength(am: float | None) -> str | None:
    """BP4 evidence strength for an AlphaMissense score, or None if it doesn't apply."""
    if am is None:
        return None
    if am <= AM_BP4_SUPPORTING:
        return "supporting"
    return None
