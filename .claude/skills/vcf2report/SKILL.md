---
name: vcf2report
description: >
  THE vcf2report command — the single guided entry point, on Claude Desktop (MCP) AND
  Claude Code / terminal. Turns an exome VCF into an auditable ACMG variant-interpretation
  report (laudo) through 8 visible stages: dependency check → inspect → annotate → analysis
  triage → prioritize (gnomAD + AlphaMissense + ClinVar + HPO) → QC → laudo (rendered inline
  as an Artifact). Use whenever the user wants to analyze a VCF / exome, get a variant report
  or "laudo", classify variants by ACMG, or install / set up vcf2report. (The acmg-classify,
  variant-report and vcf2report-orchestrator skills are lower-level Claude Desktop / MCP
  references — this skill is the whole guided flow.)
---

# vcf2report — guided ACMG analysis

You are the operator's guide. Turn one exome VCF into an **auditable ACMG laudo**, driving
**8 visible stages** from this chat. Run what you can yourself; when a step needs the user's
machine, data, or a large download, give the **exact command / tool call**, confirm, and wait
for the output before continuing. Be concise, confirm before anything networked or heavy, and
**never present the result as a diagnosis — it is a draft for expert review.**

## One engine, three renderers — detect your surface

The vcf2report engine is identical everywhere: it reads gnomAD / AlphaMissense / ClinVar / HPO
**directly on the machine** (DuckDB + tabix). Only *how you drive it* and *how you show progress*
changes. Both MCP and Bash are thin adapters over the same package — no data flows through MCP.

| | **Claude Desktop** (classic) | **Claude for Mac** | **Claude Code** / terminal |
|---|---|---|---|
| Drive the engine | vcf2report **MCP tools** | vcf2report **MCP tools** | **Bash** (`scripts/`) |
| Progress surface | plain-text ✅/⚠️ checklist | **show_widget** stepper | **Workflow → Background Tasks** pane |
| Command panel | markdown menu | show_widget | show_widget |
| Laudo | **Artifact** | **Artifact** | **Artifact** |

**Pick the path by which tools exist (best-effort):** if the **`Workflow`** tool is available →
Code path (render stages as Background-Tasks phases, drive via Bash). Else if the vcf2report
**MCP tools** (`data_status`, `inspect_vcf`, `run_report`…) are available → Desktop path (drive via
MCP, show a show_widget stepper on Mac or a text checklist on classic Desktop). If neither, run the
stages as plain text. **A stage's ACTION always runs; only its rendering degrades** — and a stage
that no-ops (already annotated, no phenotype) shows a visible "skipping" note, never a silent skip.

## Setup (once)

1. **Locate/install the project.** Look in the current dir, then `~/vcf2report`; else clone into a
   **fresh** location (⚠️ never *inside* an existing `vcf2report/` folder):
   ```bash
   git clone https://github.com/gbbarra/vcf2report.git ~/vcf2report && cd ~/vcf2report
   python3 -m pip install -e .        # only runtime dep is jinja2 — demo databases are bundled
   ```
   Confirm `pyproject.toml` + `scripts/run_headless.py` exist; dry-run `python3 scripts/run_headless.py`
   → it must print 5 candidates.
2. **Desktop only:** the vcf2report MCP server must be registered (`pip install -e ".[mcp]"`, add the
   block from `claude_desktop_config.example.json`, restart). `data_status` confirms it's live.
3. **Optional — AlphaMissense** (calibrated missense PP3/BP4, ~1 GB, confirm first):
   `bash scripts/fetch_alphamissense.sh && python3 scripts/freeze_alphamissense.py`. Without it the
   pipeline still runs (missense defers to VUS). Network is **off by default** (`VCF2REPORT_ALLOW_NETWORK=1`).

## The 8-stage flow

On **Code**, after Stages 1–2 gather the VCF, drive Stages 1&3–8 with **one Workflow call** —
`scriptPath: <repo>/.claude/skills/vcf2report/references/analyze.workflow.js`,
`args: { repo, vcf, sample, hpo, phenotypeText, out, lift, chain, reference }` — its phases render
as the Background-Tasks boxes. On **Desktop**, call the MCP tool named for each stage in order.

### 1 · 🖥️ Dependency check — FIRST; this opens the progress surface
Before asking anything else, run the readiness probe and **render it visually, explaining each item
and what it enables/disables**:
- **Code:** `python3 scripts/preflight.py` (it is also the first Workflow phase, which is what
  auto-opens Background Tasks).
- **Desktop:** the `data_status` MCP tool.

Show python ≥ 3.10, which of `bcftools`/`snpEff`/`vcfanno` are on PATH, and for **each store**
(gnomAD parquet, AlphaMissense, ClinVar, HPO) whether it is present + what it costs if missing —
**loudly** for gnomAD (missing → PM2/BA1/BS1 disabled, absence not assertable → **over-call risk**).

### 2 · Ask for the VCF + phenotype
- **VCF** — single-proband **GRCh38** `.vcf`/`.vcf.gz`. Warn on another build → set `lift: true`.
  Derive `<SAMPLE>` from the `#CHROM` header (`grep -m1 '^#CHROM' <VCF> | cut -f10-`) / `bcftools
  query -l` / filename stem.
- **Phenotype** — HPO ids, or free text (Claude maps it in Stage 6). GA4GH phenopacket →
  `python3 scripts/phenopacket_to_inputs.py <file>`.
- Render the **command panel** here: read `references/command_panel.html`, substitute `{{SAMPLE}}`
  and the build-state tokens (GRCh38 → `{{STATUS}}`=`GRCh38 · ready`,`{{STATUS_BG}}`=`bg-success`,
  `{{STATUS_FG}}`=`text-success`; GRCh37 → `GRCh37 · liftover`,`bg-warning`,`text-warning`), publish
  via `show_widget` (markdown menu if unavailable). Optionally `spawn_task` `"VCF2Report - <SAMPLE>"`
  for a per-sample session.

### 3 · 🖥️ Inspect VCF — is it annotated?
- **Code:** `python3 scripts/inspect_vcf.py <VCF> [--hpo]`  ·  **Desktop:** `inspect_vcf(vcf_path)`.

Reports build, sample, total + PASS counts, and **annotated?** (VEP CSQ / SnpEff ANN / consequence /
none) plus the capability map. This decides whether Stage 4 needs to annotate.

### 4 · 🖥️ Annotate — only if not annotated
- **Already annotated** → visible "skipping (consequence terms present)".
- **Not annotated** + `bcftools`+`snpEff`+`vcfanno` + a GRCh38 **reference FASTA** present → annotate:
  **Code** `bash scripts/annotate_vcf.sh <VCF> <REF> <OUT>.annotated.vcf.gz` · **Desktop**
  `annotate_vcf(vcf_path, reference, out_dir)`.
- **Else** explain: classification is **coordinate-only** — PVS1/PM4/PP3/BP4 and HGVS c./p. are
  unavailable; gnomAD/ClinVar coordinate lookups + the ≥2★ ClinVar safety flag still work
  (`docs/ANNOTATION.md`). Never invent HGVS.

### 5 · 🤖 Analysis triage — what can we conclude?
- **Code/Desktop:** `analysis_capabilities(vcf_path, hpo_given)` (or the `capabilities` block from
  `inspect_vcf.py`). State **each ACMG criterion → available | limited | na** with the one-line
  consequence (e.g. no gnomAD store → PM2 disabled; single-proband → segregation N/A). This is the
  **honesty gate**, stated *before* running.

### 6 · 🖥️ Prioritize — gnomAD + AlphaMissense + ClinVar + HPO
- First **map free-text phenotype → HPO** (Claude), write one `HP:xxxxxxx` per line to `<OUT>/hpo.txt`.
- Then run the engine: **Code** `python3 scripts/run_headless.py <VCF> --hpo <FILE> --sample-id
  <SAMPLE> --out <OUT> --timing` · **Desktop** `run_report(vcf_path, hpo_terms, sample_id, out_dir)`.
- Read back the **ranked candidate list** (gene → tier, phenotype/tier-topped first), the funnel
  (variants → candidates), and the per-stage **timings**. All local + deterministic — no LLM, no network.

### 7 · 🖥️ QC — the gate
Surface QC as its own step: the funnel (total → PASS → QC-passing → candidates), the
**sequencing-quality panel** (depth/GQ, Ti/Tv, het:hom, indel:SNV, multiallelic, novelty), and any
**gnomAD-store / coverage safety-net warning**. Code: from the run output / `report.qc`; Desktop:
`parse_vcf` + `run_report`'s `report.qc`.

### 8 · 🤖 Laudo — the auditable report (Artifact)
1. Read `<OUT>/<name>_report.md` (**Desktop:** use `run_report`'s inline `markdown`).
2. Build a self-contained HTML laudo and publish it with the **Artifact** tool, using
   `references/report_template.html` (in the cloned repo) as the exact style + structure. Fill:
   masthead (sample, build, generated, HPO chips), the **Conclusion**, the QC funnel, the
   **Sequencing-quality** panel, **one card per candidate** (gene, HGVS c./p. + transcript,
   coordinate, tier pill, the data-fact row, the `rule_path`, and the full ACMG criteria table with
   met / N-A / — states), then Methods + the disclaimer footer. If the VCF was not VEP/SnpEff-
   annotated, HGVS/transcript are blank (show the coordinate) — say so, don't invent HGVS.
3. Give the Artifact link + a 2–3 line summary (primary vs secondary/ACMG-SF, any ABraOM-dropped).

## Guardrails (always)
- **Compact layout.** Render short values (QC metrics, per-variant facts) **inline** (the template's
  `.kv` class) so they flow horizontally — never a tall narrow column with empty sides.
- **Draft, not diagnostic.** Keep the "not for clinical use" banner in every laudo.
- **GRCh38 only.** Flag any other build (→ liftover in Stage 3).
- **Privacy.** The VCF never leaves the machine; confirm before any network step.
- Prefer running steps yourself; hand a command to the user only when it needs their data, their
  credentials, or a large download.
