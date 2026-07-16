"""Environment readiness — the single source of truth for the dependency check.

Both the MCP ``data_status`` tool (Claude Desktop) and ``scripts/preflight.py``
(Claude Code / terminal) call :func:`readiness`, so the guided flow's Stage 1
reports identically on every surface. MCP-free by design: this module must import
without the ``mcp`` SDK so the Bash path stays lightweight.
"""
from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path

from . import config
from . import stores as _stores


def snpeff_jar() -> Path:
    """Where scripts/setup_snpeff.sh installs SnpEff (SNPEFF_JAR overrides)."""
    return Path(os.environ.get("SNPEFF_JAR") or (config.DATA_DIR / "tools" / "snpEff" / "snpEff.jar"))


def _snpeff_available() -> bool:
    return snpeff_jar().is_file() or bool(shutil.which("snpEff"))


def _store(present: bool, path, enables: str) -> dict:
    return {"present": bool(present), "path": str(path) if path else None, "enables": enables}


def readiness() -> dict:
    """What is installed and what each piece enables/disables in the ACMG run.

    Returns a superset dict: the original ``data_status`` keys (kept for
    back-compat) plus ``python``, ``stores`` (each with an ``enables`` note), so
    Stage 1 can explain to the user why a missing store matters.
    """
    # snpEff is normally a JAR (scripts/setup_snpeff.sh), not a PATH command — checking only
    # shutil.which would report it missing while annotation actually works. vcfanno is legacy:
    # the engine reads gnomAD/AlphaMissense/ClinVar from the Parquet stores, so it is not
    # required for annotation and its absence costs nothing.
    tools = {t: bool(shutil.which(t)) for t in ("bcftools", "vcfanno")}
    tools["snpEff"] = _snpeff_available()
    gnomad_parquet = config._resolve_gnomad_parquet()
    stores = {
        "gnomad_parquet": _store(
            gnomad_parquet, gnomad_parquet,
            "PM2 / BA1 / BS1 — rare/common population frequency. Missing → these are "
            "disabled and absence can't be asserted (over-call risk)."),
        "alphamissense": _store(
            config.ALPHAMISSENSE_LOCAL.exists(), config.ALPHAMISSENSE_LOCAL,
            "PP3 / BP4 — calibrated missense pathogenicity. Missing → missense defers to VUS."),
        "clinvar_tabix": _store(
            config.CLINVAR_TABIX.exists(), config.CLINVAR_TABIX,
            "PS1 / PM5 / PP5 / BP6 + the ≥2★ ClinVar safety flag. Missing → falls back to "
            "the bundled slice / live NCBI."),
        # This table now drives three things, not one: PP4, and — via the inheritance modes
        # HPO annotates per gene — the PVS1 recessive-LoF route and the PM2/BS1 AF ceilings.
        # Its absence degrades conservatively (PVS1 falls back to the constraint-only gate,
        # frequencies to the strict default) but SILENTLY, and recessive disease genes go
        # back to being unreachable by PVS1 — so name the full cost here.
        "hpo": _store(
            config.HPO_GENES_LOCAL.exists(), config.HPO_GENES_LOCAL,
            "PP4 gene↔phenotype overlap + gene mode-of-inheritance (the PVS1 recessive-LoF "
            "route and the PM2/BS1 frequency ceilings). Missing → no phenotype "
            "prioritisation, PVS1 falls back to population constraint alone (blind to "
            "recessive disease genes), and every gene takes the strict default AF ceiling."),
    }
    bundled = {
        "sample_vcf": config.SAMPLE_VCF.exists(),
        "clinvar_slice": config.CLINVAR_LOCAL.exists(),
        "gnomad_snapshot": config.GNOMAD_LOCAL.exists(),
        "abraom": config.ABRAOM_LOCAL.exists(),
        "hpo": config.HPO_GENES_LOCAL.exists(),
    }
    return {
        "python": platform.python_version(),
        "package_importable": True,  # this code is running, so the import succeeded
        "annotation_tools_on_path": tools,
        "stores": stores,
        # Per-store size + build date + freshness (quick: no row scan). A full integrity /
        # completeness scan is scripts/check_stores.py (or store_health(measure=True)).
        "store_health": _stores.store_health(measure=False),
        "bundled_local_data": bundled,
        "ready_for_offline_demo": all(bundled.values()),
        # What annotation actually requires: bcftools + snpEff. vcfanno is legacy (the engine
        # reads the Parquet stores), so gating on it would report annotation as unavailable on
        # a machine where it works fine.
        "annotation_tools_installed": tools["bcftools"] and tools["snpEff"],
        "network_egress_allowed": config.allow_network(),
        "note": "Patient data stays local by default: no gnomAD/NCBI calls unless "
                "VCF2REPORT_ALLOW_NETWORK=1. Tools on PATH do not imply the annotation "
                "databases are present; run scripts/setup_data.sh.",
    }
