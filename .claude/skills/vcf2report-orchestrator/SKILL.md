---
name: vcf2report-orchestrator
description: >-
  End-to-end clinical exome interpretation. Use when a user provides a VCF file
  (exome/genome) and wants a draft variant report, or asks to "analyze this VCF",
  "interpret these variants", or "find the candidate for this phenotype". Drives
  parsing, QC, annotation, ACMG classification, and report generation via the
  vcf2report MCP tools.
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
1. `parse_vcf(vcf_path)` → report the genome build and the QC funnel (total →
   PASS → QC-passing). If the build is not GRCh38, warn the user and stop unless
   they confirm — annotations assume GRCh38.
2. Decide candidates. The fastest path is `run_report(vcf_path, hpo_terms)`, which
   runs the whole pipeline and returns the tiers, the ABraOM-filtered list, and
   the rendered Markdown. For a guided walkthrough instead, per candidate call
   `gnomad_frequency`, `clinvar_lookup`, `abraom_frequency`, and
   `hpo_phenotype_match`, then `classify_variant`.
3. For each classified variant, invoke the **acmg-classify** skill's rules when you
   need to explain or adjudicate judgment criteria (PM1, PS3, PP2, PM5). Never
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
