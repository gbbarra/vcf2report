---
name: analyze-vcf
description: >
  Guided harness to run vcf2report on an exome VCF entirely from Claude — set up
  the environment and databases, run the pipeline, and render an auditable ACMG
  variant-interpretation report (laudo) inline as an Artifact. Use whenever the
  user wants to analyze a VCF / exome, get a variant report or "laudo", classify
  variants by ACMG, or install / set up vcf2report.
---

# vcf2report — guided analysis harness

You are the operator's guide. Turn an exome VCF into an **auditable ACMG variant
report**, driving the whole flow from this chat: run what you can yourself via
`Bash`, and when a step needs the user's machine or data, give the **exact
command**, ask them to run it in their terminal, and wait for the output before
continuing. Be concise, confirm before anything networked or heavy, and never
present the result as a diagnosis — it is a **draft for expert review**.

Work through the steps in order. Skip a step only if a check shows it is already done.

## Step 0 — Locate or install the project
This skill can be installed at the **user level** and run from anywhere, so it must
not assume the repo is already present. First put the project on the machine:
1. Look for it: the current directory, then `~/vcf2report`, then ask the user for a
   path — `test -f <dir>/pyproject.toml`.
2. If it is nowhere, clone it into a **fresh** location (⚠️ never *inside* an existing
   `vcf2report/` folder — that nests and checks out an almost-empty tree):
   ```bash
   git clone https://github.com/gbbarra/vcf2report.git ~/vcf2report && cd ~/vcf2report
   ```
3. `cd` into the repo and confirm `pyproject.toml` and `scripts/run_headless.py` exist.

## Step 1 — Environment & install
1. Check Python: `python3 --version` (need **≥ 3.10**).
2. Check the package: `python3 -c "import vcf2report; print('ok')"`.
3. If it's not installed, from the repo root:
   ```bash
   python3 -m pip install -e .          # only runtime dep is jinja2 — demo databases are bundled
   ```
   Run it yourself if you can; otherwise hand the user the command and wait.
4. Confirm with a dry run on the bundled sample: `python3 scripts/run_headless.py`
   → it must print 5 candidates. If it does, the environment is good.

## Step 2 — Optional: AlphaMissense (calibrated missense PP3/BP4)
Only if the user wants to recover pathogenic **missense** variants (v2). It needs a
~1 GB download, so **confirm first**:
```bash
brew install htslib                         # or: conda install -c bioconda htslib
bash scripts/fetch_alphamissense.sh         # download + index (~1 GB, CC BY 4.0)
python3 scripts/freeze_alphamissense.py      # local, offline
```
Without it, the pipeline still runs (missense pathogenicity defers to VUS). Live
gnomAD/ClinVar are also opt-in (`VCF2REPORT_ALLOW_NETWORK=1`) — never on by default.

## Step 3 — Collect the inputs
Ask the user for:
- **VCF** — path to a single-proband **GRCh38** `.vcf` or `.vcf.gz`. (Warn if it's
  another build; coordinate lookups are skipped then.)
- **Phenotype** — the patient's **HPO terms** (e.g. `HP:0001250 Seizure`). If they
  give free text ("seizures, developmental delay"), map it to HPO ids for them.
  Write the terms to a file, one `HP:xxxxxxx` per line:
  ```bash
  printf 'HP:0001250\nHP:0002133\n' > /tmp/hpo.txt
  ```
  (If they have a GA4GH phenopacket: `python3 scripts/phenopacket_to_inputs.py <file>`.)

## Step 3.5 — App surfaces (sample name · command panel · per-sample session)
Make the run feel like an app inside Claude Code. **Best-effort** — if a tool below is
not available in this environment, skip it silently and continue with plain text.

1. **Sample name.** Read it from the VCF — the sample columns of the `#CHROM` header
   line (`grep -m1 '^#CHROM' <VCF> | cut -f10-`), or `bcftools query -l <VCF>` if
   present; fall back to the filename stem. Use it as `<SAMPLE>` everywhere below.
2. **Per-sample session (left panel).** Offer this analysis as its own entry in
   Recents: call `spawn_task` with `title: "VCF2Report - <SAMPLE>"` and a prompt that
   re-runs `/analyze-vcf` for this VCF. One click spins it into a titled session. Skip
   if `spawn_task` isn't available.
3. **Command panel (center panel).** Render the interactive control panel: read
   `references/command_panel.html`, substitute `{{SAMPLE}}` with the sample name and
   the status tokens for the build state — GRCh38 ready → `{{STATUS}}`=`GRCh38 · pronto`,
   `{{STATUS_BG}}`=`bg-success`, `{{STATUS_FG}}`=`text-success`; GRCh37 → `GRCh37 ·
   liftover`, `bg-warning`, `text-warning` — then publish it with the `visualize`
   `show_widget` tool. Its buttons drive each step via `sendPrompt`. If `show_widget`
   isn't available, list the same actions as a short markdown menu instead.

## Step 4 — Run the analysis (visible phases, right panel)
Run the pipeline so the user watches it work in the **Background Tasks** pane. Two ways:
- **Rich (named phases)** — when there's real background work (a liftover or a gnomAD
  build): invoke the `Workflow` tool with `scriptPath:
  <repo>/.claude/skills/analyze-vcf/references/analyze.workflow.js` and `args: { repo,
  vcf, sample, hpo, out, lift, buildGnomad, jobs }`. The phases **Setup → Frequencies →
  Classify → Report** show live in the right panel. Set `lift: true` if Step 3 detected
  GRCh37; `buildGnomad: true` for offline population frequency when no local table
  exists (`build_gnomad_local.py --from-vcf --jobs 24` — an exome build is ~45–60 min,
  so confirm first and mention the ETA).
- **Lean** — just run it yourself:
  ```bash
  python3 scripts/run_headless.py <VCF> --hpo <HPO_FILE> --out <OUT_DIR> --timing
  ```
Either way, read back the **candidate list** (gene → tier), the funnel (variants →
candidates), and the per-stage **timings** in plain terms.

## Step 5 — Render the laudo (inline)
1. Read the generated report markdown (`<OUT_DIR>/<name>_report.md`).
2. Build a self-contained HTML report from it and publish it with the **Artifact**
   tool, using the template at
   `<repo>/.claude/skills/analyze-vcf/references/report_template.html` (from the repo
   located/cloned in Step 0) as the exact style + structure. If this skill is
   installed at the user level, that template still lives in the cloned repo.
   Fill it from the run: masthead (sample, build, generated, HPO chips), the
   **Conclusion** (executive summary, from the report's Conclusion section), the QC
   funnel, the **Sequencing quality** panel (depth/GQ at variant sites, Ti/Tv,
   het:hom, indel:SNV, multiallelic, novelty, PASS), one card per candidate (gene,
   **HGVS c./p. and transcript** from the report's Transcript + "Variant (c./p.)"
   columns, coordinate, tier pill, the data-fact row, the `rule_path`, and the full
   ACMG criteria table with met / N-A / — states), then Methods + the disclaimer
   footer. If the VCF was not VEP/SnpEff-annotated the HGVS/transcript will be blank
   (the report shows the coordinate instead) — say so rather than inventing HGVS.
3. Give the user the Artifact link and a 2–3 line summary of the findings
   (primary vs secondary/ACMG-SF, and any ABraOM-dropped candidates).

## Guardrails (always)
- **Compact layout.** Save vertical space: when a section holds several short values
  (QC metrics, per-variant facts), render them **inline** (the template's `.kv` class)
  so they flow horizontally and wrap to fill the width — never a tall narrow column
  with empty sides.
- **Draft, not diagnostic.** Keep the "not for clinical use" banner in every laudo.
- **GRCh38 only.** Flag any other build.
- **Privacy.** The VCF never leaves the machine; confirm before any network step.
- Prefer running steps yourself via `Bash`; hand a command to the user only when it
  needs their data, their credentials, or a large download.
