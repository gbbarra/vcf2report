"""Render a :class:`ReportModel` to Markdown (primary) or HTML.

Uses Jinja2 when available; otherwise a built-in Markdown renderer keeps the
pipeline dependency-free. Markdown renders natively inside Claude Desktop, so
the reviewer sees the draft immediately.
"""

from __future__ import annotations

from pathlib import Path

from .. import config, demo
from .assemble import ReportModel, carrier_findings, split_findings, summarize
from .vus_triage import probable_pathogenic_vus


def render_markdown(report: ReportModel) -> str:
    """Render Markdown, preferring the Jinja2 template if present."""
    template = config.TEMPLATES_DIR / "report.md.j2"
    try:
        import jinja2

        if template.exists():
            env = jinja2.Environment(
                loader=jinja2.FileSystemLoader(str(config.TEMPLATES_DIR)),
                trim_blocks=True,
                lstrip_blocks=True,
                autoescape=False,
            )
            env.filters["pct"] = lambda x: "n/a" if x is None else f"{x:.6f}"
            env.filters["kvjoin"] = kvjoin
            # Same trimming rule as the builtin path below — the two renderers must agree, or
            # whether a laudo is 6.8 MB or 1 MB depends on whether jinja2 happens to be installed.
            env.filters["informative"] = lambda crits, const=(), trim=True: [
                c
                for c in crits
                if not trim or (c.code not in const and not _uninformative(c))
            ]
            env.tests["uninformative"] = _uninformative
            primary, secondary, other = split_findings(report.classifications)
            carriers = carrier_findings(report.classifications)
            # Carriers are routed into `other` by split_findings; pull them out so the
            # template can give them their own section instead of burying a reproductively
            # relevant finding among unrelated VUS.
            other = [c for c in other if c not in carriers]
            vus = probable_pathogenic_vus(report.classifications)
            return env.get_template("report.md.j2").render(
                r=report,
                primary=primary,
                secondary=secondary,
                other=other,
                carriers=carriers,
                vus=vus,
                conclusion=summarize(report),
                run_constant=run_constant_criteria(report.classifications),
                trim=len(report.classifications) >= HOIST_MIN_VARIANTS,
            )
    except ImportError:
        pass
    return _render_markdown_builtin(report)


def _uninformative(cr) -> bool:
    """A criterion row that neither fired nor recorded anything: no evidence, no citation.

    Its Reasoning is then a generic statement of why the criterion could not apply here — real
    information, but the SAME information for every variant it is printed against. Contrast a
    row that did not fire but DID record evidence (PM2 with ``gnomad_af=0.0714``): that one
    proves the engine looked and what it saw, which is exactly what an audit needs.
    """
    return (not cr.met) and (not cr.evidence) and (not cr.citation)


#: Below this many classifications, no criterion is hoisted — see run_constant_criteria.
HOIST_MIN_VARIANTS = 10


def run_constant_criteria(classifications) -> dict:
    """Criteria that fired for NO variant in the run, mapped to a summary of why.

    A criterion that never fires across an entire exome is a fact about the ANALYSIS, not about
    any allele. A single-proband VCF cannot supply PS2 / PM3 / PM6 / PP1 / BP2 / BS4 whatever the
    variant is; PS3 / BS3 / BP3 / BP5 / PS4 wait on evidence this run did not gather; and
    PS1 / PM1 / PM5 sit idle whenever the ClinVar residue index was never built. Repeating that
    once per variant does not make it more auditable — it makes it unfindable.

    Measured on SYN-016 (1,178 classifications, 6.76 MB): 10 criteria produced byte-identical
    rows 1,178 times each, 1.42 MB / 21% of the report. Seven more never fired either but carried
    a per-variant evidence dict — PM1 alone cost 356 KB printing `window=7, cutoff=3,
    enrichment_cutoff=2.0` for 1,178 variants it never assessed, because the residue index was
    absent. Those are parameters, not observations.

    The reason is usually identical everywhere; where it is not (PM1 has 12 variants of it,
    mostly "not a missense variant"), the dominant one is named with its count and the remainder
    is acknowledged rather than flattened. Nothing is dropped from ``_results.json``.
    """
    # The hoist exists to remove REPETITION. Below a handful of variants there is none to
    # remove, and "fired for no variant" becomes trivially true of almost every criterion —
    # on a single-variant report it would empty the table entirely and call that a summary.
    # Whole rationale section at this size is tens of KB; print it in full.
    if len(classifications) < HOIST_MIN_VARIANTS:
        return {}
    from collections import Counter

    reasons, order = {}, []
    for c in classifications:
        for cr in c.criteria:
            if cr.code not in reasons:
                reasons[cr.code] = Counter()
                order.append(cr.code)
            if cr.met:
                reasons[cr.code] = None  # fired at least once -> not run-constant
            elif reasons[cr.code] is not None:
                reasons[cr.code][cr.reasoning] += 1

    out = {}
    for code in order:
        counts = reasons[code]
        if counts is None or not counts:
            continue
        (top, n), rest = counts.most_common(1)[0], len(counts) - 1
        out[code] = (
            top
            if rest == 0
            else (
                f"{top} — for {n} of {sum(counts.values())} variants; "
                f"{rest} other reason(s) apply to the remainder (see the JSON)"
            )
        )
    return out


def kvjoin(evidence) -> str:
    """Render a criterion's evidence dict for the human-readable trail.

    A missing value prints as an em dash, never as ``None``. Python's ``None`` is a single
    token standing for at least six different things across this report — the index was not
    built, the criterion does not apply to this consequence, the source was not consulted, the
    upstream value was itself unavailable, the record was found but carries no such field, and
    the metric does not exist for this gene. Measured on SYN-073: 2,202 evidence fields, all
    printing the same word.

    The Reasoning column already distinguishes all six correctly — the engine knows which is
    which. Evidence's job is the concrete VALUE, and the truthful value here is "there was
    none". Printing a Python literal in a clinical document also invites a reader to parse it
    as data. The JSON keeps `null`, which is what a machine should see.
    """
    if not evidence:
        return "—"
    return ", ".join(f"{k}={'—' if v is None else v}" for k, v in evidence.items())


def _zyg(v) -> str:
    """Zygosity for display. A half-call is shown as an assumption, not an observation."""
    z = v.zygosity or "?"
    return f"{z} (partial call)" if getattr(v, "partial_call", False) else z


def _fmt_af(x) -> str:
    return "n/a" if x is None else f"{x:.6f}"


def _render_markdown_builtin(report: ReportModel) -> str:
    L: list[str] = []
    L.append(f"# Variant Interpretation Report — {report.sample_id}")
    L.append("")
    L.append(
        "> **DRAFT — for expert review. Not for clinical use.** "
        "Auto-generated candidate interpretation to be verified and signed out "
        "by a qualified professional."
    )
    L.append("")
    # The demo stamp is repeated here as a masthead banner, not only inside the conclusion
    # bullets, because this renderer is what MCP and --stdout emit and a reader may see only
    # the first screen.
    if report.provenance.get("mode") == "demo":
        L.append("> 🧪 " + demo.stamp_line(report.provenance))
        L.append("")
    L.append(f"- **Genome build:** {report.build}")
    L.append(f"- **Pipeline:** vcf2report v{report.tool_version}")
    L.append(f"- **Generated:** {report.generated}")
    L.append(
        f"- **Patient HPO terms:** {', '.join(report.hpo_terms) or 'none provided'}"
    )
    L.append("")

    L.append("## Conclusion (draft interpretation)")
    L.append("")
    for line in summarize(report):
        L.append(f"- {line}")
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
    if q.near_splice_excluded:
        L.append(
            f"- Set aside by design: **{q.near_splice_excluded}** rare splice-*adjacent* "
            "variants (3–8 bases into the intron, or the last exonic bases). The canonical "
            "±1,2 splice sites ARE included; these are not, because no splice predictor is "
            "wired in. A splice-region variant already known to ClinVar as P/LP is still "
            "shortlisted. **This is a stated sensitivity limit of the shortlist.**"
        )
    for w in q.warnings:
        L.append(f"- ⚠️ {w}")
    L.append("")
    if q.qc_rescued:
        # These variants are NOT in any candidate table — the QC gate removed them before
        # annotation. Without this section the report simply would not contain them.
        L.append("### Removed by QC, but known to ClinVar")
        L.append("")
        L.append(
            "Variants the per-variant QC gate dropped for a borderline metric that "
            "ClinVar nevertheless classifies Pathogenic / Likely Pathogenic. They were "
            "**not classified** and are **not** candidates below — they are listed so a "
            "marginal-quality call on a known allele is confirmed, not silently lost:"
        )
        for note in q.qc_rescued:
            L.append(f"- ⚠️ {note}")
        L.append("")
    if q.local_cohort_filtered:
        L.append("### Local-cohort frequency filtering")
        L.append("")
        L.append(
            "Spurious candidates a gnomAD-only pipeline would have kept, removed "
            "using the operator's local cohort frequencies:"
        )
        for note in q.local_cohort_filtered:
            L.append(f"- {note}")
        L.append("")

    sq = report.seq_quality
    if sq:
        L.append("## Sequencing quality (estimated from variant sites)")
        L.append("")
        L.append(
            f"- **Assay (by variant count):** {sq.assay_guess} "
            f"({sq.n_variants} variants)"
        )
        if sq.dp_mean is not None:
            L.append(
                f"- **Depth at called sites:** {sq.dp_mean}x mean / "
                f"{sq.dp_median}x median — {sq.dp_pct_ge10}% ≥10x, "
                f"{sq.dp_pct_ge20}% ≥20x"
            )
        if sq.gq_median is not None:
            L.append(
                f"- **Genotype quality:** median {sq.gq_median}, {sq.gq_pct_ge20}% ≥20"
            )
        if sq.titv is not None:
            L.append(f"- **Ti/Tv (SNVs):** {sq.titv} ({sq.n_snv} SNVs)")
        if sq.het_hom_ratio is not None:
            L.append(
                f"- **Het/Hom:** {sq.het_hom_ratio} ({sq.n_het} het / {sq.n_hom} hom)"
            )
        if sq.indel_snv_ratio is not None:
            L.append(
                f"- **Indel:SNV ratio:** {sq.indel_snv_ratio} ({sq.n_indel} indels)"
            )
        if sq.pct_multiallelic is not None:
            L.append(
                f"- **Multiallelic sites:** {sq.pct_multiallelic}% "
                f"({sq.n_multiallelic_sites} / {sq.n_sites})"
            )
        if sq.pct_novel is not None:
            L.append(
                f"- **Novel (not in dbSNP):** {sq.pct_novel}% "
                f"({sq.n_with_rsid} carry an rsID)"
            )
        if sq.pct_het_ab_balanced is not None:
            L.append(
                f"- **Het allele balance:** {sq.pct_het_ab_balanced}% balanced "
                f"({config.QC_AB_MIN:g}–{config.QC_AB_MAX:g})"
            )
        if sq.pct_pass is not None:
            L.append(f"- **FILTER = PASS:** {sq.pct_pass}%")
        for note in sq.notes:
            L.append(f"- _{note}_")
        L.append("")

    primary, secondary, other = split_findings(report.classifications)

    def _findings_table(rows):
        if not rows:
            L.append("_None._")
            L.append("")
            return
        L.append(
            "| Gene | Transcript | Variant (c./p.) | Zyg | Consequence | ClinVar | gnomAD AF | local cohort AF | HPO | ACMG |"
        )
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for c in rows:
            v, a = c.variant, c.annotation
            hgvs = " ".join(x for x in [v.hgvs_c, v.hgvs_p] if x) or v.key
            L.append(
                f"| {v.gene or '?'} | {v.transcript or '—'} | {hgvs} | {_zyg(v)} | {v.consequence or '?'} "
                f"| {a.clinvar_significance or '—'} | {_fmt_af(a.gnomad_af)} | {_fmt_af(a.local_cohort_af)} "
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
    L.append(
        "_P/LP variants in ACMG SF v3.2 genes, unrelated to the indication — "
        "reportable actionable secondary findings, subject to the patient's "
        "opt-in policy._"
    )
    L.append("")
    _findings_table(secondary)

    carriers = carrier_findings(report.classifications)
    if carriers:
        L.append("## Carrier findings (recessive, heterozygous)")
        L.append("")
        L.append(
            "_Pathogenic/Likely-Pathogenic alleles found in the heterozygous state in genes "
            "whose only known disease mechanism is **recessive**. The ACMG tier refers to the "
            "VARIANT, not to this patient: a single copy does **not** cause the condition and "
            "does **not** explain the indication — the patient is a healthy carrier. Most "
            "people carry a few. Listed here for **reproductive relevance** (partner testing / "
            "genetic counselling), deliberately kept out of the diagnostic sections so they "
            "cannot compete with the indication._"
        )
        L.append("")
        _findings_table(carriers)

    other_non_carrier = [c for c in other if c not in carriers]
    if other_non_carrier:
        L.append("## Other candidates")
        L.append("")
        L.append(
            "_Incidental P/LP not on the ACMG SF list, plus phenotype-unrelated "
            "uncertain/benign candidates. Not routinely reported._"
        )
        L.append("")
        _findings_table(other_non_carrier)

    vus = probable_pathogenic_vus(report.classifications)
    if vus:
        L.append("## Probable-pathogenic VUS — for expert review")
        L.append("")
        L.append(
            "_Phenotype-relevant variants the engine held at **Uncertain Significance** (correctly — "
            "the evidence is Supporting-only) but which carry suggestive molecular signal. The ACMG "
            "tier is unchanged; these are ranked for a human + model second look (literature on the "
            "residue/gene, domain / functional context, splicing prediction, and the ClinVar "
            "assertion's underlying evidence). Claude can help work through each._"
        )
        L.append("")
        L.append("| Rank | Gene | Variant (c./p.) | Suggestive evidence |")
        L.append("|---|---|---|---|")
        for i, e in enumerate(vus, 1):
            v = e["classification"].variant
            hgvs = " ".join(x for x in [v.hgvs_c, v.hgvs_p] if x) or v.key
            sig = "; ".join(
                f"{s['signal']} ({s['value']})"
                if s["value"] is not True
                else s["signal"]
                for s in e["signals"]
            )
            L.append(f"| {i} | {v.gene or '?'} | {hgvs} | {sig} |")
        L.append("")

    L.append("## Per-variant ACMG rationale (auditable)")
    # Hoist what is true of the whole run, and say what the tables leave out. Auditability is
    # the point of this section, so the omission is stated, never silent — and the JSON keeps
    # all 28 criteria per variant regardless of what is rendered here.
    constant = run_constant_criteria(report.classifications)
    trim = len(report.classifications) >= HOIST_MIN_VARIANTS
    if constant:
        L.append("")
        L.append(
            f"{len(constant)} criteria fired for **no** variant in this analysis — a fact "
            "about the analysis rather than about any allele. They are stated once here and "
            "omitted from the tables below:"
        )
        L.append("")
        for code, why in constant.items():
            L.append(f"- **{code}** — {why}")
        L.append("")
        L.append(
            "Tables also omit any further criterion that neither fired nor recorded "
            "evidence. A criterion that did not fire but *did* record what it saw is kept, "
            "because that is the line proving the engine looked. The complete 28-criterion "
            "set for every variant is in the accompanying `_results.json`."
        )
    for c in report.classifications:
        v = c.variant
        L.append("")
        tx = f"{v.transcript}:" if v.transcript else ""
        hp = f" ({v.hgvs_p})" if v.hgvs_p else ""
        label = f"{tx}{v.hgvs_c}{hp}" if (v.hgvs_c or v.hgvs_p) else v.key
        L.append(f"### {v.gene or '?'} — {label} → {c.tier}")
        L.append("")
        L.append(f"**Rule path:** `{c.rule_path}`")
        L.append("")
        L.append(
            "| Criterion | Applied | Strength | Evidence | Source | By | Reasoning |"
        )
        L.append("|---|---|---|---|---|---|---|")
        omitted = 0
        for cr in c.criteria:
            if trim and cr.code in constant:
                continue  # covered by the run-wide statement above
            if trim and _uninformative(cr):
                omitted += 1  # counted per block, never silently dropped
                continue
            if not cr.applies:
                state = "N/A"
            else:
                state = "✅ met" if cr.met else "—"
            strength = cr.applied_strength or cr.default_strength
            evidence = kvjoin(cr.evidence)  # same rule as the Jinja path — never "None"
            source = "; ".join(cr.citation) or "—"
            L.append(
                f"| **{cr.code}** | {state} | {strength} | {evidence} | {source} "
                f"| {cr.adjudicated_by} | {cr.reasoning} |"
            )
        if omitted:
            L.append("")
            L.append(
                f"_{omitted} further criterion(s) neither fired nor recorded evidence for "
                "this variant; see the JSON for the full set._"
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
    L.append(
        "- Single-proband analysis: criteria requiring parental/segregation/"
        "phasing data (PS2, PM3, PM6, PP1, BP2, BS4) are reported as N/A."
    )
    L.append(
        "- Judgment criteria (PS3, PS4, BS3, BP3, BP5) are surfaced for expert/model "
        "adjudication and default to not-met unless explicitly supported."
    )
    L.append(
        "- Population and clinical databases are versioned snapshots; re-check "
        "before sign-out."
    )
    L.append("- **This is a draft-generation aid, not a diagnostic device.**")
    L.append("")
    return "\n".join(L)


def write_report(report: ReportModel, out_dir: Path | None = None) -> Path:
    """Write the Markdown laudo and its companion ``<sample>_results.json``.

    The JSON is not optional output: it is the queryable form of the same run (every variant
    with its full ACMG criterion trail, the routed buckets, the conclusion, the ClinVar
    do-not-dismiss list), and it is what makes a follow-up conversation a lookup instead of a
    re-analysis. Writing it here rather than at each call site is what keeps the two in step —
    only the CLI used to write it, so the MCP path produced a laudo that could not be explored
    at all, which is the whole point of persisting the run.

    Returns the Markdown path; the JSON sits beside it.
    """
    out_dir = out_dir or config.OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    md = render_markdown(report)
    fp = out_dir / f"{report.sample_id}_report.md"
    fp.write_text(md)
    from .explore import write_explore  # local import: avoids a cycle via assemble

    write_explore(report, out_dir / f"{report.sample_id}_results.json")
    return fp


def results_json_for(report_path: Path | str) -> Path:
    """The ``_results.json`` companion of a ``_report.md`` path."""
    return Path(str(report_path).replace("_report.md", "_results.json"))
