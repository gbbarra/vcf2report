"""End-to-end pipeline: VCF path + HPO terms -> auditable ReportModel.

This is the single orchestration entry point shared by the headless CLI and the
MCP ``run_report`` tool. Each stage is thin and delegates to the module that
owns it (parse / qc / annotate / filter / acmg / report).
"""
from __future__ import annotations

from pathlib import Path

from . import config
from .acmg.engine import classify
from .annotate import annotate_variant
from .models import Classification, QCSummary
from .report.assemble import ReportModel, build_report
from .vcf.filter import filter_variants
from .vcf.parse import parse_vcf
from .vcf.qc import apply_qc


def run_pipeline(
    vcf_path: str | Path,
    hpo_terms: list[str] | None = None,
    sample_id: str | None = None,
    max_af: float = config.AF_RECESSIVE_MAX,
) -> ReportModel:
    hpo_terms = hpo_terms or []
    vcf_path = Path(vcf_path)
    sample_id = sample_id or vcf_path.stem.replace(".vcf", "")

    variants, build, _header = parse_vcf(vcf_path)
    qc = QCSummary(total_variants=len(variants), build=build or "unknown")

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

    annotated = [(v, annotate_variant(v, hpo_terms, build_trusted=build_trusted))
                 for v in kept]
    candidates, funnel = filter_variants(annotated, max_af=max_af)
    qc.after_rarity = funnel.after_rarity
    qc.after_impact = funnel.after_impact
    qc.candidates = funnel.candidates
    qc.abraom_filtered = funnel.abraom_filtered

    classifications: list[Classification] = [classify(v, a) for v, a in candidates]
    return build_report(sample_id, hpo_terms, qc, classifications)
