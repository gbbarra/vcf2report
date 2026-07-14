"""Inspect a VCF and derive which analyses it supports — Stages 3 & 5 of the flow.

Pure functions shared by the MCP tools (``inspect_vcf`` / ``analysis_capabilities``,
Claude Desktop) and the terminal path (``scripts/inspect_vcf.py``, Claude Code):
detect whether the VCF is annotated, sniff the build and sample, and map that plus
the installed stores to the ACMG criteria that are actually computable — so the run
is honest about what it can and cannot conclude.
"""
from __future__ import annotations

import gzip
from pathlib import Path

from . import config
from .status import readiness
from .vcf.parse import parse_vcf


def _open_text(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _annotation_info_keys() -> set:
    keys = {"ANN", "CSQ"}
    for aliases in config.INFO_ALIASES.values():
        keys.update(aliases)
    return keys


def _sample_id_from_header(vcf_path: str) -> str | None:
    try:
        with _open_text(vcf_path) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    cols = line.rstrip("\n").split("\t")
                    return cols[9] if len(cols) > 9 else None
                if not line.startswith("#"):
                    break
    except OSError:
        return None
    return None


def _detect_annotation(variants):
    """Return (annotated: bool, source: str|None). Consequence field wins; then the
    SnpEff ANN / VEP CSQ / population INFO blocks. Samples the first 200 records."""
    sample = variants[:200]
    if any(v.consequence for v in sample):
        return True, "consequence"
    info_keys = _annotation_info_keys()
    for name, present in (("SnpEff ANN", ("ANN",)), ("VEP CSQ", ("CSQ",))):
        if any(any(k in (v.info or {}) for k in present) for v in sample):
            return True, name
    if any(any(k in (v.info or {}) for k in info_keys) for v in sample):
        return True, "population INFO"
    return False, None


def inspect_vcf(vcf_path: str) -> dict:
    """Build, sample, variant counts, and annotation status for one VCF (Stage 3)."""
    variants, build, _ = parse_vcf(vcf_path)
    annotated, source = _detect_annotation(variants)
    pass_filter = sum(1 for v in variants if v.filter_status in ("PASS", ".", "", None))
    present_keys = sorted({k for v in variants[:200] for k in (v.info or {})
                           if k in _annotation_info_keys()})
    return {
        "vcf_path": str(vcf_path),
        "build": build,
        "sample_id": _sample_id_from_header(vcf_path) or Path(vcf_path).stem,
        "total_variants": len(variants),
        "pass_filter": pass_filter,
        "annotated": annotated,
        "annotation_source": source,
        "info_keys_present": present_keys,
    }


_AVAILABLE, _LIMITED, _NA = "available", "limited", "na"


def analysis_capabilities(vcf_path: str, hpo_given: bool = False,
                          inspection: dict | None = None, rd: dict | None = None) -> dict:
    """Map the VCF + installed stores to each ACMG criterion's status (Stage 5).

    Each entry is {status: available|limited|na, reason}. ``available`` = the data
    to evaluate it is present; ``limited`` = evaluable but degraded/absent-data;
    ``na`` = not applicable to this input (e.g. segregation on a single proband).
    """
    insp = inspection or inspect_vcf(vcf_path)
    rd = rd or readiness()
    ann = insp.get("annotated")
    build_ok = insp.get("build") == config.GENOME_BUILD if hasattr(config, "GENOME_BUILD") else True
    st = rd.get("stores", {})
    has_gnomad = st.get("gnomad_parquet", {}).get("present")
    has_am = st.get("alphamissense", {}).get("present")
    has_clinvar = st.get("clinvar_tabix", {}).get("present") or rd.get(
        "bundled_local_data", {}).get("clinvar_slice")
    has_hpo = st.get("hpo", {}).get("present")

    def crit(status, reason):
        return {"status": status, "reason": reason}

    caps = {}
    ann_reason = ("consequence annotation present" if ann
                  else "VCF not annotated (no VEP/SnpEff/consequence) — annotate first (Stage 4)")
    caps["PVS1 (LoF)"] = crit(_AVAILABLE if ann else _LIMITED, ann_reason)
    caps["PM4 (in-frame / stop-loss)"] = crit(_AVAILABLE if ann else _LIMITED, ann_reason)
    caps["PP3 / BP4 (missense)"] = crit(
        _AVAILABLE if (ann and has_am) else _LIMITED,
        "AlphaMissense present" if has_am else "AlphaMissense store absent → missense defers to VUS")
    caps["PM2 / BA1 / BS1 (frequency)"] = crit(
        _AVAILABLE if has_gnomad else _LIMITED,
        "local gnomAD store detected" if has_gnomad
        else "gnomAD store absent → frequency criteria disabled, absence not assertable (over-call risk)")
    caps["PS1 / PM5 / PP5 / BP6 (ClinVar)"] = crit(
        _AVAILABLE if has_clinvar else _LIMITED,
        "ClinVar store present" if has_clinvar else "no local ClinVar → live NCBI only (network-gated)")
    caps["PP4 (phenotype)"] = crit(
        _AVAILABLE if (has_hpo and hpo_given) else (_LIMITED if has_hpo else _NA),
        "HPO terms supplied" if (has_hpo and hpo_given)
        else "no phenotype given — genotype-only run" if has_hpo
        else "HPO store absent")
    caps["PS2 / PM3 / PM6 / PP1 / BS4 (segregation)"] = crit(
        _NA, "single-proband input — de novo / in-trans / segregation need trio or family data")
    if not build_ok:
        caps["_build_warning"] = crit(
            _LIMITED, f"build is {insp.get('build')}, not {getattr(config, 'GENOME_BUILD', 'GRCh38')} "
                      "— coordinate lookups skipped until lifted over")
    return {"annotated": ann, "annotation_source": insp.get("annotation_source"),
            "build": insp.get("build"), "criteria": caps}
