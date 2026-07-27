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
from .annotate import add_alphamissense, add_clinvar_residue, annotate_variant
from .models import Classification, QCSummary
from .report.assemble import ReportModel, build_report
from .vcf import seqqc
from .vcf.filter import filter_variants
from .vcf.parse import parse_vcf
from .vcf.qc import apply_qc


# A QC rescue re-annotates variants the gate dropped, so it is bounded: a real exome
# yields ~100 of them, but a pathological input must not turn into an unbounded second
# annotation pass. If the cap bites, the report says so rather than implying full coverage.
QC_RESCUE_MAX = 500


def _qc_loss_warnings(qc: QCSummary, variants: list, kept: list,
                      dropped: list[tuple]) -> None:
    """Warn when the QC gate, or the annotation upstream of it, silently ate the callset.

    QC is the one stage that can remove EVERYTHING and still produce a clean-looking
    report: the funnel prints ``0 candidates`` and the conclusion states there is no
    finding. Each guard below turns a specific silent-loss mode into a loud warning, in
    the same spirit as the gnomAD-parquet guard further down.
    """
    if not variants:
        qc.warnings.append(
            "The VCF contained no variant records — nothing was analysed. Check that the "
            "file is a variant callset and not an empty or header-only VCF."
        )
        return

    reasons: dict[str, int] = {}
    for _v, reason in dropped:
        key = reason.split("=")[0].split(" (")[0]
        reasons[key] = reasons.get(key, 0) + 1
    non_carrier = reasons.get("non-carrier", 0)

    if not kept:
        detail = ", ".join(f"{k}: {n}" for k, n in sorted(reasons.items(), key=lambda x: -x[1]))
        qc.warnings.append(
            f"ALL {len(variants)} variants were removed by per-variant QC ({detail}) — the "
            "report below has nothing to analyse and its 'no finding' statements carry NO "
            "evidential weight. Check the genotype/quality fields before reading further."
        )
    if non_carrier and non_carrier >= 0.9 * len(variants):
        qc.warnings.append(
            f"{non_carrier} of {len(variants)} variants ({non_carrier / len(variants):.0%}) were "
            "dropped as non-carriers. That is the signature of a sites-only VCF (no FORMAT/sample "
            "column) or of genotypes written as './.' — this pipeline analyses ONE proband and "
            "needs that proband's genotypes. Re-export the VCF with the sample column."
        )
    if kept and not any(v.consequence for v in kept):
        qc.warnings.append(
            f"None of the {len(kept)} post-QC variants carry a consequence annotation "
            "(SnpEff ANN / VEP CSQ / a plain CSQ key), so the impact filter can only drop "
            "them: the candidate list will be empty or near-empty for a reason that is NOT "
            "biological. Annotate the VCF (SnpEff or VEP) before running."
        )


def _qc_rescue(qc: QCSummary, dropped: list[tuple], hpo_terms: list[str],
               build_trusted: bool) -> None:
    """Name QC-dropped variants that ClinVar classifies P/LP with criteria-backed review.

    The report's do-not-dismiss net (``assemble.clinvar_pathogenic_flags``) only sees
    variants that survived to classification, and QC runs before annotation — so a
    well-reviewed pathogenic allele one GQ point under threshold was deleted from the
    report entirely while the conclusion announced 'no finding'. Only borderline
    quantitative drops are reconsidered (see ``qc.is_metric_drop``); the ACMG tier is
    NOT computed and nothing is re-admitted to the candidate list. This is a flag.

    The bar is >=1 star (ClinVar "criteria provided"), one step below the report's
    do-not-dismiss net. The net triages among variants the reader can already SEE in the
    ranked table; a QC drop is invisible, so a criteria-backed pathogenic assertion is
    worth naming even when only a single submitter made it.
    """
    from .report.assemble import clinvar_stars
    from .vcf.filter import is_impactful
    from .vcf.qc import is_metric_drop

    pool = [(v, r) for v, r in dropped if is_metric_drop(r) and is_impactful(v.consequence)]
    if not pool:
        return
    capped = len(pool) > QC_RESCUE_MAX
    for v, reason in pool[:QC_RESCUE_MAX]:
        # Resolve ClinVar through the SAME path the classified variants use (VCF INFO
        # first, then the local/live client). Calling the DB client directly would miss
        # every pre-annotated exome — the recommended production input.
        a = annotate_variant(v, [], build_trusted=build_trusted,
                             with_alphamissense=False, with_clinvar_residue=False)
        sig = (a.clinvar_significance or "").lower().replace("_", " ")
        stars = clinvar_stars(a.clinvar_review_status)
        if not (sig.startswith("pathogenic") or sig.startswith("likely pathogenic")):
            continue
        if stars < 1:
            continue
        qc.qc_rescued.append(
            f"{v.gene or v.key} {v.hgvs_p or v.hgvs_c or v.key} — ClinVar "
            f"{a.clinvar_significance} ({stars}★) but dropped at QC ({reason}). NOT classified "
            "and NOT counted as a candidate; confirm the call orthogonally (Sanger / "
            "re-sequencing) before dismissing it."
        )
    if capped:
        qc.warnings.append(
            f"The QC-drop rescue examined the first {QC_RESCUE_MAX} of {len(pool)} borderline "
            f"coding/splice drops; {len(pool) - QC_RESCUE_MAX} were NOT checked against ClinVar."
        )


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

    kept, dropped = apply_qc(variants)
    qc.after_qc = len(kept)
    _qc_loss_warnings(qc, variants, kept, dropped)
    _mark("qc_s")

    # gnomAD frequency is needed for the rarity filter across the WHOLE post-QC set, so
    # if a DuckDB/Parquet store is configured, resolve them all in one vectorised join
    # up front — annotate_variant's per-variant gnomad.lookup then reads that cache
    # instead of ~11k tabix/remote round-trips. No-op when the parquet isn't configured.
    from .annotate import gnomad_parquet
    primed = gnomad_parquet.prime(kept)
    # Safety net: if the operator pointed us at a parquet store but 0 of a real
    # callset resolved, the store is unavailable (drive unmounted) or empty. Left
    # silent, every variant looks absent from gnomAD: the rarity filter can no longer
    # exclude common variants (unknown AF is treated as 0) and BA1/BS1 can't fire to
    # down-weight them -> gross over-calling. Flag it loudly, don't ship a wrong report.
    if config.GNOMAD_PARQUET and len(kept) >= 50 and primed == 0:
        from .annotate.gnomad_parquet import _get_duckdb
        store_present = Path(config.GNOMAD_PARQUET).exists()
        if store_present and _get_duckdb() is None:
            cause = "the 'duckdb' package is not installed — run `pip install duckdb`"
        elif not store_present:
            cause = ("the store is unavailable (e.g. the drive is unmounted) or empty — "
                     "mount it and re-run")
        else:
            cause = "the store matched no variants (a schema or coordinate mismatch)"
        qc.warnings.append(
            f"gnomAD parquet is configured but 0 of {len(kept)} post-QC variants "
            f"resolved from it — {cause}. Population frequencies were NOT applied: the "
            "rarity filter cannot exclude common variants and BA1/BS1 cannot down-weight "
            "them, so the report likely OVER-calls."
        )
    _mark("gnomad_prime_s")
    # ClinVar — the one source with no batch path — resolved in ONE chr-pruned DuckDB join
    # over the post-QC set when a Parquet store is present; clinvar.lookup reads that cache
    # first, then falls back to the per-variant tabix / live / slice. No-op if unconfigured.
    from .annotate import clinvar_parquet
    from .vcf.filter import is_impactful as _impactful
    from .vcf.qc import is_metric_drop as _metric_drop
    # The QC-drop rescue queries ClinVar for variants that never reach annotation, so they
    # must ride along in the same primed join rather than falling back to per-variant lookups.
    _rescue_pool = [v for v, r in dropped if _metric_drop(r) and _impactful(v.consequence)]
    clinvar_parquet.prime(kept + _rescue_pool[:QC_RESCUE_MAX])
    _qc_rescue(qc, dropped, hpo_terms, build_trusted)
    _mark("clinvar_prime_s")
    # AlphaMissense is deferred: it only feeds PP3/BP4 at classification, never the
    # filter, so we skip the (per-variant, ~1 GB tabix) lookup across the whole
    # post-QC set and query just the surviving candidates below.
    annotated = [(v, annotate_variant(v, hpo_terms, build_trusted=build_trusted,
                                      with_alphamissense=False, with_clinvar_residue=False))
                 for v in kept]
    _mark("annotate_s")
    candidates, funnel = filter_variants(annotated, max_af=max_af)
    qc.after_rarity = funnel.after_rarity
    qc.after_impact = funnel.after_impact
    qc.candidates = funnel.candidates
    qc.abraom_filtered = funnel.abraom_filtered
    qc.near_splice_excluded = funnel.near_splice_excluded
    _mark("filter_s")

    if build_trusted:
        from .annotate import alphamissense
        alphamissense.prime([v for v, a in candidates])
        for v, a in candidates:
            add_alphamissense(v, a)
    _mark("alphamissense_s")

    # Residue-level ClinVar evidence (PS1/PM5/PM1) — deferred here for the same reason as
    # AlphaMissense: it feeds classification only, never the filter, so the surviving candidates
    # are the only variants that need it.
    for v, a in candidates:
        add_clinvar_residue(v, a, build_trusted=build_trusted)
    _mark("clinvar_residue_s")

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
