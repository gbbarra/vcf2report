"""Render a :class:`ReportModel` to Markdown (primary) or HTML.

Uses Jinja2 when available; otherwise a built-in Markdown renderer keeps the
pipeline dependency-free. Markdown renders natively inside Claude Desktop, so
the reviewer sees the draft immediately.
"""
from __future__ import annotations

from pathlib import Path

from .. import config
from .assemble import ReportModel, split_findings


def render_markdown(report: ReportModel) -> str:
    """Render Markdown, preferring the Jinja2 template if present."""
    template = config.TEMPLATES_DIR / "report.md.j2"
    try:
        import jinja2

        if template.exists():
            env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(str(config.TEMPLATES_DIR)),
                trim_blocks=True, lstrip_blocks=True, autoescape=False,
            )
            env.filters["pct"] = lambda x: "n/a" if x is None else f"{x:.6f}"
            env.filters["kvjoin"] = lambda d: (
                ", ".join(f"{k}={v}" for k, v in d.items()) or "—"
            )
            primary, secondary, other = split_findings(report.classifications)
            return env.get_template("report.md.j2").render(
                r=report, primary=primary, secondary=secondary, other=other)
    except ImportError:
        pass
    return _render_markdown_builtin(report)


def _fmt_af(x) -> str:
    return "n/a" if x is None else f"{x:.6f}"


def _render_markdown_builtin(report: ReportModel) -> str:
    L: list[str] = []
    L.append(f"# Variant Interpretation Report — {report.sample_id}")
    L.append("")
    L.append("> **DRAFT — for expert review. Not for clinical use.** "
             "Auto-generated candidate interpretation to be verified and signed out "
             "by a qualified professional.")
    L.append("")
    L.append(f"- **Genome build:** {report.build}")
    L.append(f"- **Pipeline:** vcf2report v{report.tool_version}")
    L.append(f"- **Generated:** {report.generated}")
    L.append(f"- **Patient HPO terms:** {', '.join(report.hpo_terms) or 'none provided'}")
    L.append("")

    q = report.qc
    L.append("## Quality control & filtering funnel")
    L.append("")
    L.append(f"- Total variants: **{q.total_variants}**")
    L.append(f"- PASS filter: **{q.pass_filter}**")
    L.append(f"- After QC (DP/GQ/AB): **{q.after_qc}**")
    L.append(f"- After rarity: **{q.after_rarity}**")
    L.append(f"- After coding/splice impact: **{q.after_impact}**")
    L.append(f"- **Candidates classified: {q.candidates}**")
    for w in q.warnings:
        L.append(f"- ⚠️ {w}")
    L.append("")
    if q.abraom_filtered:
        L.append("### Brazilian-frequency filtering (ABraOM)")
        L.append("")
        L.append("Spurious candidates a gnomAD-only pipeline would have kept, removed "
                 "using ABraOM (SABE) local frequencies:")
        for note in q.abraom_filtered:
            L.append(f"- {note}")
        L.append("")

    sq = report.seq_quality
    if sq:
        L.append("## Sequencing quality (estimated from variant sites)")
        L.append("")
        L.append(f"- **Assay (by variant count):** {sq.assay_guess} "
                 f"({sq.n_variants} variants)")
        if sq.dp_mean is not None:
            L.append(f"- **Depth at called sites:** {sq.dp_mean}x mean / "
                     f"{sq.dp_median}x median — {sq.dp_pct_ge10}% ≥10x, "
                     f"{sq.dp_pct_ge20}% ≥20x")
        if sq.gq_median is not None:
            L.append(f"- **Genotype quality:** median {sq.gq_median}, "
                     f"{sq.gq_pct_ge20}% ≥20")
        if sq.titv is not None:
            L.append(f"- **Ti/Tv (SNVs):** {sq.titv} ({sq.n_snv} SNVs)")
        if sq.het_hom_ratio is not None:
            L.append(f"- **Het/Hom:** {sq.het_hom_ratio} "
                     f"({sq.n_het} het / {sq.n_hom} hom)")
        for note in sq.notes:
            L.append(f"- _{note}_")
        L.append("")

    primary, secondary, other = split_findings(report.classifications)

    def _findings_table(rows):
        if not rows:
            L.append("_None._")
            L.append("")
            return
        L.append("| Gene | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | ABraOM AF | HPO | ACMG |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for c in rows:
            v, a = c.variant, c.annotation
            hgvs = " ".join(x for x in [v.hgvs_c, v.hgvs_p] if x) or v.key
            L.append(
                f"| {v.gene or '?'} | {hgvs} | {v.zygosity or '?'} | {v.consequence or '?'} "
                f"| {a.clinvar_significance or '—'} | {_fmt_af(a.gnomad_af)} | {_fmt_af(a.abraom_af)} "
                f"| {a.hpo_match_score if a.hpo_match_score is not None else '—'} | **{c.tier}** |"
            )
        L.append("")

    L.append("## Primary (diagnostic) findings")
    L.append("")
    L.append("_Variants in genes overlapping the patient's phenotype._")
    L.append("")
    _findings_table(primary)

    L.append("## Secondary findings (ACMG SF v3.2)")
    L.append("")
    L.append("_P/LP variants in ACMG SF v3.2 genes, unrelated to the indication — "
             "reportable actionable secondary findings, subject to the patient's "
             "opt-in policy._")
    L.append("")
    _findings_table(secondary)

    if other:
        L.append("## Other candidates")
        L.append("")
        L.append("_Incidental P/LP not on the ACMG SF list, plus phenotype-unrelated "
                 "uncertain/benign candidates. Not routinely reported._")
        L.append("")
        _findings_table(other)

    L.append("## Per-variant ACMG rationale (auditable)")
    for c in report.classifications:
        v = c.variant
        L.append("")
        L.append(f"### {v.gene or '?'} — {v.hgvs_p or v.hgvs_c or v.key} → {c.tier}")
        L.append("")
        L.append(f"**Rule path:** `{c.rule_path}`")
        L.append("")
        L.append("| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |")
        L.append("|---|---|---|---|---|---|---|")
        for cr in c.criteria:
            if not cr.applies:
                state = "N/A"
            else:
                state = "✅ met" if cr.met else "—"
            strength = cr.applied_strength or cr.default_strength
            evidence = ", ".join(f"{k}={v2}" for k, v2 in cr.evidence.items()) or "—"
            source = "; ".join(cr.citation) or "—"
            L.append(
                f"| **{cr.code}** | {state} | {strength} | {evidence} | {source} "
                f"| {cr.adjudicated_by} | {cr.reasoning} |"
            )
    L.append("")

    L.append("## Methods")
    L.append("")
    for k, val in report.methods.items():
        L.append(f"- **{k}:** {val}")
    L.append("")

    if report.timings:
        L.append("## Performance (this run)")
        L.append("")
        for k, val in report.timings.items():
            unit = "" if k in ("variants_per_s",) else " s"
            L.append(f"- **{k.replace('_s', '').replace('_', ' ')}:** {val}{unit}")
        L.append("")

    L.append("## Limitations & disclaimers")
    L.append("")
    L.append("- Single-proband analysis: criteria requiring parental/segregation/"
             "phasing data (PS2, PM3, PM6, PP1, BS4) are reported as N/A.")
    L.append("- Judgment criteria (PS3, PS4, PM1, PM5, PP2) are surfaced for expert/model "
             "adjudication and default to not-met unless explicitly supported.")
    L.append("- Population and clinical databases are versioned snapshots; re-check "
             "before sign-out.")
    L.append("- **This is a draft-generation aid, not a diagnostic device.**")
    L.append("")
    return "\n".join(L)


def write_report(report: ReportModel, out_dir: Path | None = None) -> Path:
    out_dir = out_dir or config.OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    md = render_markdown(report)
    fp = out_dir / f"{report.sample_id}_report.md"
    fp.write_text(md)
    return fp
