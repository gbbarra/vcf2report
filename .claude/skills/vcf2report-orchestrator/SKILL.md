---
name: vcf2report-orchestrator
description: >-
  Claude Desktop only — REQUIRES the vcf2report MCP server. End-to-end clinical
  exome interpretation via the vcf2report MCP tools (parse, QC, annotate, ACMG,
  report). Do NOT use in plain Claude Code (the MCP tools are absent there) — the
  `analyze-vcf` skill runs the whole pipeline via the terminal instead.
---

# vcf2report — clinical exome interpretation SOP

You are assisting a laboratory professional. Turn a raw exome VCF into a **draft**
ACMG variant report. The output is always a draft for expert review — never a
final clinical result.

## Inputs to confirm first
1. **VCF path** — ask for it if not given.
2. **Patient phenotype** — HPO terms (e.g. `HP:0001250`) or free-text you map to
   HPO. Phenotype drives candidate ranking and the PP4 criterion, so do not skip it.

## Steps (use the MCP tools of the `vcf2report` server)
0. On a new machine, call `data_status()` to confirm what's ready (bundled data
   runs the demo; the annotation tools are needed for a raw real exome).
1. `parse_vcf(vcf_path)` → report the genome build and the QC funnel (total →
   PASS → QC-passing). If the build is not GRCh38, warn the user and stop unless
   they confirm — annotations assume GRCh38.
2. Produce the report. The one-call path is `annotate_and_report(vcf_path,
   hpo_terms, reference)`: it annotates a raw VCF locally (SnpEff + vcfanno) if
   needed, else classifies an already-annotated VCF directly, and returns tiers,
   the ABraOM-filtered list, per-stage timings, and the Markdown. (If the VCF is
   already annotated you can also call `run_report(vcf_path, hpo_terms)`.) For a
   guided walkthrough of one candidate, call `gnomad_frequency`, `clinvar_lookup`,
   `abraom_frequency`, `hpo_phenotype_match`, then `classify_variant`.
3. For each classified variant, invoke the **acmg-classify** skill's rules when you
   need to explain or adjudicate judgment criteria (PS3, PS4, PM1, PM5, PP2). Never
   invent evidence; every criterion you assert must cite a source + accession + date.
4. Invoke the **variant-report** skill to present the draft, keeping the auditable
   ACMG criterion tables verbatim.
5. Surface the ABraOM-filtered variants explicitly — they are candidates a
   gnomAD-only pipeline would have wrongly kept.

## Guardrails
- Label every output **"DRAFT — for expert review, not for clinical use."**
- Do not upgrade a tier beyond what the combining rules return.
- State N/A criteria honestly (single-proband analysis cannot assess trio /
  segregation / phasing criteria).
- If a tool returns no data for a variant, say so rather than guessing.
