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

import shutil
import subprocess
from pathlib import Path
from typing import Optional

from . import config
from .acmg.engine import classify as _classify
from .annotate import annotate_variant
from .annotate import abraom as _abraom
from .annotate import clinvar as _clinvar
from .annotate import gnomad as _gnomad
from .annotate import hpo as _hpo
from .models import Variant
from .pipeline import run_pipeline
from .report.render import render_markdown, write_report
from .vcf.parse import parse_vcf as _parse_vcf
from .vcf.qc import apply_qc

try:
    from mcp.server.fastmcp import FastMCP
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "The MCP server needs the MCP SDK. Install it with: pip install 'mcp[cli]'"
    ) from exc

mcp = FastMCP("vcf2report")


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
def abraom_frequency(chrom: str, pos: int, ref: str, alt: str) -> dict:
    """ABraOM (Brazilian SABE cohort) allele frequency — local population check."""
    return _abraom.lookup(_mk_variant(chrom, pos, ref, alt))


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
    """
    report = run_pipeline(vcf_path, hpo_terms=hpo_terms or [], sample_id=sample_id or None)
    out = Path(out_dir) if out_dir else config.OUTPUT_DIR
    fp = write_report(report, out)
    return {
        "report_path": str(fp),
        "candidates": report.qc.candidates,
        "abraom_filtered": report.qc.abraom_filtered,
        "tiers": [{"gene": c.variant.gene, "variant": c.variant.hgvs_p or c.variant.key,
                   "tier": c.tier, "rule_path": c.rule_path} for c in report.classifications],
        "markdown": render_markdown(report),
    }


@mcp.tool()
def data_status() -> dict:
    """Report readiness: annotation tools on PATH + bundled local datasets.

    Call this first on a new machine to see what's ready. The bundled data runs
    the offline demo immediately; the annotation tools + downloaded databases are
    needed to annotate a raw real exome (see docs/SETUP.md, docs/ANNOTATION.md).
    """
    tools = {t: bool(shutil.which(t)) for t in ("bcftools", "snpEff", "vcfanno")}
    bundled = {
        "sample_vcf": config.SAMPLE_VCF.exists(),
        "clinvar_slice": config.CLINVAR_LOCAL.exists(),
        "gnomad_snapshot": config.GNOMAD_LOCAL.exists(),
        "abraom": config.ABRAOM_LOCAL.exists(),
        "hpo": config.HPO_GENES_LOCAL.exists(),
    }
    return {
        "annotation_tools_on_path": tools,
        "bundled_local_data": bundled,
        "ready_for_offline_demo": all(bundled.values()),
        # Tools present ≠ ready to annotate: the databases (gnomAD/ClinVar VCFs,
        # SnpEff DB, reference FASTA) must also be downloaded — see docs/SETUP.md.
        "annotation_tools_installed": all(tools.values()),
        # Network egress is OFF by default; live gnomAD/ClinVar lookups only happen
        # when VCF2REPORT_ALLOW_NETWORK=1 (and OFFLINE is not set).
        "network_egress_allowed": config.allow_network(),
        "note": "Patient data stays local by default: no gnomAD/NCBI calls unless "
                "VCF2REPORT_ALLOW_NETWORK=1. Tools on PATH do not imply the annotation "
                "databases are present; run scripts/setup_data.sh.",
    }


def _annotation_info_keys() -> set:
    # Keep in sync with the reader: every INFO alias + the consequence blocks.
    keys = {"ANN", "CSQ"}
    for aliases in config.INFO_ALIASES.values():
        keys.update(aliases)
    return keys


def _looks_annotated(variants) -> bool:
    sample = variants[:200]
    if any(v.consequence for v in sample):
        return True
    keys = _annotation_info_keys()
    return any(any(k in (v.info or {}) for k in keys) for v in sample)


@mcp.tool()
def annotate_and_report(vcf_path: str, hpo_terms: Optional[list[str]] = None,
                        reference: str = "", out_dir: str = "") -> dict:
    """One-call path for a bench scientist: raw VCF -> draft report.

    If the VCF is already annotated (ANN/CSQ or population fields in INFO) it is
    classified directly. Otherwise, if the annotation tools and a GRCh38
    ``reference`` FASTA are available, it is annotated locally (SnpEff + vcfanno)
    first. Returns the report path, tiers, timings, and the rendered Markdown.
    """
    variants, _build, _ = _parse_vcf(vcf_path)
    used = vcf_path
    steps: list[str] = []
    if not _looks_annotated(variants):
        have_tools = all(shutil.which(t) for t in ("bcftools", "snpEff", "vcfanno"))
        if reference and have_tools:
            out = Path(out_dir) if out_dir else config.OUTPUT_DIR
            out.mkdir(parents=True, exist_ok=True)
            annotated = out / (Path(vcf_path).stem.replace(".vcf", "") + ".annotated.vcf.gz")
            script = str(config.REPO_ROOT / "scripts" / "annotate_vcf.sh")
            proc = subprocess.run(
                ["bash", script, vcf_path, reference, str(annotated)],
                capture_output=True, text=True)
            if proc.returncode != 0:
                tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-15:]
                return {"error": "annotation_failed", "exit_code": proc.returncode,
                        "stderr_tail": "\n".join(tail),
                        "hint": "See docs/ANNOTATION.md; check tool versions and the "
                                "GRCh38 reference / vcfanno.conf.toml data paths."}
            used = str(annotated)
            steps.append("annotated locally via bcftools norm + SnpEff + vcfanno")
        else:
            steps.append(
                "VCF is not annotated and (tools or reference) are unavailable; "
                "classified on coordinates only — consequence-based criteria are "
                "limited. See docs/ANNOTATION.md to annotate first.")
    result = run_report(used, hpo_terms=hpo_terms, out_dir=out_dir)
    result["steps"] = steps
    result["annotated_input"] = used
    return result


def main() -> None:  # pragma: no cover
    mcp.run()


if __name__ == "__main__":  # pragma: no cover
    main()
