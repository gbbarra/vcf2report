"""MCP server exposing the vcf2report pipeline to Claude Desktop.

Thin adapter: every tool delegates to the importable ``vcf2report`` package and
returns compact JSON-friendly dicts — big per-variant data stays on disk and we
pass keys/paths. The one deliberate exception is ``run_report``, which also
returns the finished Markdown report inline so the reviewer sees it immediately
in Claude Desktop (the report is the product). Run with::

    python -m vcf2report.mcp_server

and register it in claude_desktop_config.json (see the example in the repo root).
Requires ``pip install "mcp[cli]"``.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from . import config
from . import inspect as _inspect
from . import status as _status
from .acmg.engine import classify as _classify
from .annotate import annotate_variant
from .annotate import local_cohort as _local_cohort
from .annotate import clinvar as _clinvar
from .annotate import gnomad as _gnomad
from .annotate import hpo as _hpo
from .models import Variant
from .pipeline import run_pipeline
from .report.render import render_markdown, results_json_for, write_report
from .vcf.parse import parse_vcf as _parse_vcf
from .vcf.qc import apply_qc

# The SDK renamed the server class in v2: `mcp.server.fastmcp.FastMCP` was removed in favour
# of `mcp.server.MCPServer`. The decorator API (@mcp.tool()) is unchanged, so accept either
# and let the installed SDK decide. Without this the module raised SystemExit on any machine
# that resolved a v2 release — including CI, which does not pin.
try:
    from mcp.server.fastmcp import FastMCP as _Server        # SDK v1
except ImportError:
    try:
        from mcp.server import MCPServer as _Server          # SDK v2
    except ImportError as exc:  # pragma: no cover
        raise SystemExit(
            "The MCP server needs the MCP SDK. Install it with: pip install 'mcp[cli]'"
        ) from exc

mcp = _Server("vcf2report")


def _mk_variant(chrom: str, pos: int, ref: str, alt: str, gene: Optional[str] = None,
                consequence: Optional[str] = None, hgvs_c: Optional[str] = None,
                hgvs_p: Optional[str] = None, zygosity: Optional[str] = None) -> Variant:
    return Variant(chrom=chrom, pos=pos, ref=ref, alt=alt, gene=gene,
                   consequence=consequence, hgvs_c=hgvs_c, hgvs_p=hgvs_p,
                   zygosity=zygosity)


@mcp.tool()
def parse_vcf(vcf_path: str) -> dict:
    """Parse and QC a VCF. Returns build, counts, and the QC-passing variants.

    Use this first to confirm the input looks right before annotating.
    """
    variants, build, _ = _parse_vcf(vcf_path)
    kept, dropped = apply_qc(variants)
    pass_filter = sum(1 for v in variants if v.filter_status in ("PASS", ".", "", None))
    return {
        "build": build,
        "total_variants": len(variants),
        "pass_filter": pass_filter,
        "qc_passing": len(kept),
        "qc_dropped": [{"key": v.key, "reason": r} for v, r in dropped],
        "variants": [
            {"key": v.key, "gene": v.gene, "consequence": v.consequence,
             "hgvs_c": v.hgvs_c, "hgvs_p": v.hgvs_p, "zygosity": v.zygosity}
            for v in kept
        ],
    }


@mcp.tool()
def gnomad_frequency(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """gnomAD popmax allele frequency, allele counts, and homozygote count."""
    return _gnomad.lookup(_mk_variant(chrom, pos, ref, alt))


@mcp.tool()
def clinvar_lookup(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """ClinVar clinical significance, review status, accession, and condition."""
    return _clinvar.lookup(_mk_variant(chrom, pos, ref, alt))


@mcp.tool()
def local_cohort_frequency(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """Allele frequency from the operator-supplied local cohort, if one is configured."""
    return _local_cohort.lookup(_mk_variant(chrom, pos, ref, alt))


@mcp.tool()
def hpo_phenotype_match(gene: str, hpo_terms: list[str]) -> dict:
    """Overlap between the patient's HPO terms and a gene's known phenotypes."""
    return _hpo.match(gene, hpo_terms)


@mcp.tool()
def classify_variant(
    chrom: str, pos: int, ref: str, alt: str,
    gene: str = "", consequence: str = "", hgvs_c: str = "", hgvs_p: str = "",
    zygosity: str = "", hpo_terms: Optional[list[str]] = None,
) -> dict:
    """Annotate one variant and run the auditable ACMG classification.

    Returns the 5-tier call, the rule path, and every criterion with its
    evidence, source, and reasoning — the auditable trail for sign-out.
    """
    v = _mk_variant(chrom, pos, ref, alt, gene or None, consequence or None,
                    hgvs_c or None, hgvs_p or None, zygosity or None)
    ann = annotate_variant(v, hpo_terms or [])
    return _classify(v, ann).to_dict()


@mcp.tool()
def run_report(vcf_path: str, hpo_terms: Optional[list[str]] = None,
               sample_id: str = "", out_dir: str = "") -> dict:
    """Run the full pipeline (parse -> QC -> annotate -> filter -> ACMG -> report).

    Writes a Markdown draft report and returns its path, the tier summary, and
    the rendered Markdown for immediate review. This is the end-to-end tool.

    Also writes `results_json` — the same run in queryable form. Pass that path to the
    explore_* tools to answer follow-up questions ("why did X get PM2?", "what would it
    take to call this?") without re-running the pipeline.
    """
    report = run_pipeline(vcf_path, hpo_terms=hpo_terms or [], sample_id=sample_id or None)
    out = Path(out_dir) if out_dir else config.OUTPUT_DIR
    fp = write_report(report, out)
    return {
        "report_path": str(fp),
        "results_json": str(results_json_for(fp)),
        "candidates": report.qc.candidates,
        "local_cohort_filtered": report.qc.local_cohort_filtered,
        "tiers": [{"gene": c.variant.gene, "variant": c.variant.hgvs_p or c.variant.key,
                   "tier": c.tier, "rule_path": c.rule_path} for c in report.classifications],
        "markdown": render_markdown(report),
    }


@mcp.tool()
def data_status() -> dict:
    """Report readiness: annotation tools on PATH + local stores (Stage 1).

    Call this first on a new machine to see what's ready and what each store
    enables/disables in the ACMG run. The bundled data runs the offline demo
    immediately; the annotation tools + downloaded databases are needed to
    annotate a raw real exome (see docs/SETUP.md, docs/ANNOTATION.md).
    """
    return _status.readiness()


@mcp.tool()
def check_stores(vcf_path: str | None = None, demo: bool = False) -> dict:
    """Full health + integrity scan of the annotation Parquet stores (gnomAD, AlphaMissense,
    ClinVar): presence, size, row count, integrity (reads cleanly), completeness vs the build
    manifest, build date + source version, and freshness by cadence (ClinVar weekly; gnomAD v4.1
    / AlphaMissense frozen). ``data_status`` carries a quick summary; this does the row scan.

    Pass ``demo=True`` with ``vcf_path`` pointing at one of the repository's committed
    ``data/example/`` fixtures to get the gate's demo verdict — the exemption that lets the guided
    flow be demonstrated without the full stores. It is REFUSED (``mode="refused"``, ``ready``
    False) for any other path, and the resulting laudo is stamped either way."""
    from . import demo as _demo_mod
    from . import stores as _stores_mod
    if not (vcf_path or demo):
        return _stores_mod.store_health(measure=True)
    return _demo_mod.gate(vcf_path, demo_requested=demo, measure=True)


@mcp.tool()
def inspect_vcf(vcf_path: str) -> dict:
    """Detect build, sample, variant counts, and whether the VCF is annotated (Stage 3).

    Returns annotated (bool) + annotation_source (VEP CSQ / SnpEff ANN / consequence /
    population INFO / null) so the flow knows whether to annotate before classifying.
    """
    return _inspect.inspect_vcf(vcf_path)


@mcp.tool()
def analysis_capabilities(vcf_path: str, hpo_given: bool = False) -> dict:
    """Which ACMG criteria are computable for this VCF (Stage 5 — the honest gate).

    Combines the VCF's annotation status with the installed stores and returns each
    criterion as available | limited | na with the reason (e.g. gnomAD store absent ->
    PM2/BA1/BS1 limited; single-proband -> segregation N/A).
    """
    return _inspect.analysis_capabilities(vcf_path, hpo_given=hpo_given)


def _annotate_vcf(vcf_path: str, reference: str = "", out_dir: str = "") -> dict:
    """Detect-then-annotate (shared by the annotate_vcf tool and the express path)."""
    variants, _build, _ = _parse_vcf(vcf_path)
    annotated, source = _inspect._detect_annotation(variants)
    if annotated:
        return {"annotated_path": vcf_path, "already_annotated": True,
                "annotation_source": source, "validated": True, "missing": [],
                "steps": [f"already annotated ({source}) — skipping annotation"]}
    # vcfanno is NOT required: frequencies/clinical data come from the Parquet stores, not
    # from a vcfanno pass. A GRCh38 reference is optional too — it only enables indel
    # left-alignment. SnpEff may be a JAR (scripts/setup_snpeff.sh) rather than on PATH.
    snpeff_jar = Path(os.environ.get("SNPEFF_JAR") or (config.DATA_DIR / "tools" / "snpEff" / "snpEff.jar"))
    has_snpeff = snpeff_jar.is_file() or bool(shutil.which("snpEff"))
    missing_tools = ([] if shutil.which("bcftools") else ["bcftools (not on PATH)"]) + \
                    ([] if has_snpeff else ["snpEff (run: bash scripts/setup_snpeff.sh)"])
    if not missing_tools:
        out = Path(out_dir) if out_dir else config.OUTPUT_DIR
        out.mkdir(parents=True, exist_ok=True)
        annotated_out = out / (Path(vcf_path).stem.replace(".vcf", "") + ".annotated.vcf.gz")
        script = str(config.REPO_ROOT / "scripts" / "annotate_vcf.sh")
        cmd = ["bash", script, vcf_path, str(annotated_out)] + ([reference] if reference else [])
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        except subprocess.TimeoutExpired:
            return {"error": "annotation_timeout",
                    "hint": "Annotation exceeded 1h; check input size / tool health."}
        if proc.returncode != 0:
            tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-15:]
            return {"error": "annotation_failed", "exit_code": proc.returncode,
                    "stderr_tail": "\n".join(tail),
                    "hint": "See docs/ANNOTATION.md; check the SnpEff database and that the "
                            "VCF's chromosome naming was resolved."}
        return {"annotated_path": str(annotated_out), "already_annotated": False,
                "validated": True, "missing": [],
                "steps": ["annotated locally via bcftools norm + SnpEff (MANE)"]}
    return {"annotated_path": vcf_path, "already_annotated": False,
            "validated": False, "missing": missing_tools,
            "steps": ["not annotated and (tools or reference) unavailable; classification "
                      "will be coordinate-only — PVS1/PM4/PP3/BP4 and HGVS are limited."],
            "hint": "See docs/ANNOTATION.md to annotate first."}


@mcp.tool()
def annotate_vcf(vcf_path: str, reference: str = "", out_dir: str = "") -> dict:
    """Annotate a raw VCF locally when possible, else report the coordinate-only limits (Stage 4).

    If already annotated -> returns it unchanged (skipping). If not, and bcftools + snpEff
    are present -> annotates via scripts/annotate_vcf.sh (gene + consequence + HGVS on MANE
    transcripts). ``reference`` is optional and only enables indel left-alignment.
    Otherwise returns validated=false with what's missing.
    """
    return _annotate_vcf(vcf_path, reference=reference, out_dir=out_dir)


@mcp.tool()
def annotate_and_report(vcf_path: str, hpo_terms: Optional[list[str]] = None,
                        reference: str = "", out_dir: str = "") -> dict:
    """One-call express path for a bench scientist: raw VCF -> draft report.

    Annotates first if needed (see annotate_vcf), then runs the full pipeline.
    The guided /vcf2report flow calls the stages separately; this is the fast path.
    """
    ann = _annotate_vcf(vcf_path, reference=reference, out_dir=out_dir)
    if ann.get("error"):
        return ann
    used = ann["annotated_path"]
    result = run_report(used, hpo_terms=hpo_terms, out_dir=out_dir)
    result["steps"] = ann["steps"]
    result["annotated_input"] = used
    return result


# ---------------------------------------------------------------------------
# Explore — question-answering over a finished run, with no re-analysis
#
# `run_report` writes `<sample>_results.json` beside the laudo: the whole run in queryable
# form. These tools read it back, so a follow-up question ("why did BRCA1 get PM2?", "what
# would it take to call this?") costs a file read rather than a second pipeline pass — and
# the answers are the SAME values the laudo was rendered from, not a re-derivation that
# could disagree with it.
# ---------------------------------------------------------------------------
def _explore_load(results_json: str):
    from .report.explore import load_explore
    return load_explore(results_json)


@mcp.tool()
def explore_case(results_json: str) -> dict:
    """Open a finished run: sample, build, candidate count, the routed buckets, the
    interpretive conclusion, and the ClinVar do-not-dismiss safety list.

    Start here after run_report. `results_json` is the `<sample>_results.json` written
    beside the Markdown laudo.
    """
    from .report.explore import overview
    return overview(_explore_load(results_json))


@mcp.tool()
def explore_gene(results_json: str, gene: str, criterion: str = "") -> dict:
    """Everything the run concluded about one gene — or WHY it got one ACMG criterion.

    Without `criterion`: the gene digest (each classified variant's tier, rule path, met
    criteria, HGVS, key annotations, and which bucket it routed to).

    With `criterion` (e.g. "PM2", "PS1"): that criterion's full trail per variant — whether
    it applied and was met, at what strength, the concrete evidence values, the citation,
    the one-line reasoning, and whether the ENGINE decided it or it is left for model/expert
    adjudication. This is the audit answer, verbatim from the persisted run.
    """
    from .report.explore import criterion_basis, explain
    data = _explore_load(results_json)
    if criterion:
        return {"gene": gene, "criterion": criterion.upper(),
                "basis": criterion_basis(data, gene, criterion)}
    return explain(data, gene)


@mcp.tool()
def explore_missing_evidence(results_json: str, gene: str = "") -> dict:
    """"What would it take to call this?" — the question a conservative engine invites.

    Re-runs the published ACMG combining rules with ONE hypothetical extra line of evidence
    at each strength and reports the WEAKEST addition that would change the tier, in either
    direction, then names the criteria that could supply it and what each concretely
    requires — a trio, a functional assay, segregation data.

    Reach for this on a VUS, which is where the engine leaves most real cases. Nothing is
    re-classified: it reads the persisted trail and reports what the combining rules would
    do. Omit `gene` for every classified variant.
    """
    from .report.explore import missing_evidence
    return {"cases": missing_evidence(_explore_load(results_json), gene or None)}


@mcp.tool()
def explore_evidence_sources(results_json: str, view: str = "clinvar",
                             include_not_met: bool = False) -> dict:
    """Which findings rest on which class of evidence — for weighing circularity.

    `view="clinvar"`: findings with a met criterion drawing on a ClinVar assertion, either
    the variant's own record (PP5/BP6) or a residue-level cross-match to a different record
    (PS1/PM5). Use this to see how much of a call is ClinVar repeating itself.

    `view="missense"`: findings carrying gene- or residue-level missense evidence
    (PP2/BP1 from gnomAD missense constraint; PS1/PM5/PM1 from the ClinVar residue index) —
    the criteria that decide most missense calls. Set `include_not_met` for the honest audit
    view, which also shows why they did NOT fire, including when the residue index was never
    built for that gene.
    """
    from .report.explore import findings_citing_clinvar, missense_evidence
    data = _explore_load(results_json)
    v = (view or "clinvar").lower()
    if v == "missense":
        return {"view": v, "findings": missense_evidence(data, met_only=not include_not_met)}
    if v != "clinvar":
        return {"error": f"unknown view {view!r}; expected 'clinvar' or 'missense'"}
    return {"view": v, "findings": findings_citing_clinvar(data)}


def main() -> None:  # pragma: no cover
    mcp.run()


if __name__ == "__main__":  # pragma: no cover
    main()
