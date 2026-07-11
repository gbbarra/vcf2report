"""End-to-end pipeline: VCF path + HPO terms -> auditable ReportModel.

This is the single orchestration entry point shared by the headless CLI and the
MCP ``run_report`` tool. Each stage is thin and delegates to the module that
owns it (parse / qc / annotate / filter / acmg / report).
"""
from __future__ import annotations

import time
from pathlib import Path

from . import config
from .acmg.engine import classify
from .annotate import add_alphamissense, annotate_variant
from .models import Classification, QCSummary
from .report.assemble import ReportModel, build_report
from .vcf import seqqc
from .vcf.filter import filter_variants
from .vcf.parse import parse_vcf
from .vcf.qc import apply_qc


def run_pipeline(
    vcf_path: str | Path,
    hpo_terms: list[str] | None = None,
    sample_id: str | None = None,
    max_af: float = config.AF_RECESSIVE_MAX,
    sample: str | None = None,
) -> ReportModel:
    hpo_terms = hpo_terms or []
    vcf_path = Path(vcf_path)
    sample_id = sample_id or vcf_path.stem.replace(".vcf", "")

    timings: dict[str, float] = {}
    _t = time.perf_counter()

    def _mark(stage: str) -> None:
        nonlocal _t
        now = time.perf_counter()
        timings[stage] = round(now - _t, 4)
        _t = now

    variants, build, header = parse_vcf(vcf_path, sample=sample)
    _mark("parse_s")
    qc = QCSummary(total_variants=len(variants), build=build or "unknown")

    # Multi-sample guard: we analyse ONE proband. Warn loudly if a multi-sample
    # VCF was passed without naming the proband (we default to the first column).
    from .vcf.parse import _sample_names
    names = _sample_names(header)
    if len(names) > 1 and sample is None:
        qc.warnings.append(
            f"Multi-sample VCF ({len(names)} samples: {', '.join(names)}); analysed "
            f"the FIRST column ({names[0]}). Pass the proband's sample name to be sure."
        )

    # Genome-build guard: everything downstream assumes GRCh38. A *confirmed*
    # different build is not trusted for coordinate-keyed annotation (skip the
    # GRCh38 DBs); an undeclared build is assumed GRCh38 with a warning.
    build_trusted = True
    if build and build != config.GENOME_BUILD:
        build_trusted = False
        qc.warnings.append(
            f"VCF build detected as {build}, expected {config.GENOME_BUILD}; "
            "coordinate-based annotation (gnomAD/ClinVar/ABraOM) was SKIPPED — "
            "re-lift to GRCh38 before clinical use."
        )
    if build is None:
        qc.warnings.append(
            f"Genome build not declared in header; assuming {config.GENOME_BUILD}."
        )

    qc.pass_filter = sum(
        1 for v in variants if v.filter_status in ("PASS", ".", "", None)
    )

    kept, _dropped = apply_qc(variants)
    qc.after_qc = len(kept)
    _mark("qc_s")

    # gnomAD frequency is needed for the rarity filter across the WHOLE post-QC set, so
    # if a DuckDB/Parquet store is configured, resolve them all in one vectorised join
    # up front — annotate_variant's per-variant gnomad.lookup then reads that cache
    # instead of ~11k tabix/remote round-trips. No-op when the parquet isn't configured.
    from .annotate import gnomad_parquet
    gnomad_parquet.prime(kept)
    # AlphaMissense is deferred: it only feeds PP3/BP4 at classification, never the
    # filter, so we skip the (per-variant, ~1 GB tabix) lookup across the whole
    # post-QC set and query just the surviving candidates below.
    annotated = [(v, annotate_variant(v, hpo_terms, build_trusted=build_trusted,
                                      with_alphamissense=False))
                 for v in kept]
    _mark("annotate_s")
    candidates, funnel = filter_variants(annotated, max_af=max_af)
    qc.after_rarity = funnel.after_rarity
    qc.after_impact = funnel.after_impact
    qc.candidates = funnel.candidates
    qc.abraom_filtered = funnel.abraom_filtered
    _mark("filter_s")

    if build_trusted:
        for v, a in candidates:
            add_alphamissense(v, a)
    _mark("alphamissense_s")

    classifications: list[Classification] = [classify(v, a) for v, a in candidates]
    _mark("classify_s")

    report = build_report(sample_id, hpo_terms, qc, classifications)
    # Sequencing-quality estimate over ALL called variants (pre-filter callset).
    report.seq_quality = seqqc.estimate(variants)
    total = round(sum(timings.values()), 4)
    timings["total_s"] = total
    if total > 0:
        timings["variants_per_s"] = round(len(variants) / total, 1)
    report.timings = timings
    return report
