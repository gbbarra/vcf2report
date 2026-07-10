---
name: variant-report
description: >-
  Assemble and present a clinical variant interpretation report from classified
  variants. Use when producing the final draft report, formatting reportable
  findings, or writing per-variant ACMG rationale into the lab's template.
---

# Variant report assembly

Produce the draft report from classified variants. Prefer the engine's renderer
(`run_report` MCP tool writes Markdown to disk and returns it) so the layout and
disclaimers stay consistent; adjust wording to the lab's style in
`references/report_style.md`.

## Must-haves
- **Header:** de-identified sample ID, genome build, pipeline + database
  versions/dates, and a prominent **"DRAFT — not for clinical use"** banner.
- **Clinical context:** indication + patient HPO terms.
- **QC & filtering funnel:** the variant counts at each step, including the
  **ABraOM (Brazilian-frequency) filtering** callout — name the variants a
  gnomAD-only pipeline would have wrongly kept.
- **Reportable findings table:** gene, HGVS (c. and p.), zygosity, consequence,
  ClinVar, gnomAD AF, ABraOM AF, phenotype match, final ACMG tier.
- **Per-variant ACMG rationale:** include the auditable criterion table
  **verbatim** — code, met/not/N-A, strength, evidence value, source, engine-vs-
  model, reasoning — plus the combining-rule path to the tier. This is the point
  of the tool; never summarise it away.
- **Methods + Limitations & disclaimers:** thresholds, database versions, the
  single-proband N/A caveats, and the research-use statement.

## Tone
Concise and clinical. Do not overstate certainty. Flag where expert review is
most needed (VUS, judgment criteria, conflicting evidence).
