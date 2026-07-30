---
description: Guided ACMG exome analysis — one VCF to an auditable laudo in 8 visible stages
argument-hint: "[VCF path] [HPO file or HP: ids or free-text phenotype]"
---

Invoke the **`vcf2report` skill** now — call the `Skill` tool with `skill: "vcf2report"` — and
drive its 8-stage flow. The skill at `.claude/skills/vcf2report/SKILL.md` is the single source of
truth for the flow; do not restate or improvise it here.

The user typed: `$ARGUMENTS`

Read whatever is there as the VCF path and/or the phenotype (an HPO file, a comma-separated list of
`HP:` ids, or free text to map in Stage 6). Anything missing, ask for in Stage 2 — do not guess a
path and do not substitute a bundled sample for a file the user meant to supply.

Two things that decide whether this run is legitimate, both already specified in the skill — honour
them exactly:

- **The Stage-1 Parquet store gate is hard.** If gnomAD / AlphaMissense / ClinVar are not all
  available and intact, STOP. Do not analyse a real exome with blind stores.
- **The one exemption is demo mode**, and only for a VCF committed under `data/example/`. Pass
  `demo: true` to the analysis Workflow for those, and never for anything else — the gate refuses
  it anyway. The laudo is stamped DEMONSTRATION RUN regardless; never drop that banner.

Never present the result as a diagnosis. It is a draft for expert review.
