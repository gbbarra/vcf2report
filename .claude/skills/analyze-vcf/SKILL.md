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
   python -m pip install -e .          # only runtime dep is jinja2 — demo databases are bundled
   ```
   Run it yourself if you can; otherwise hand the user the command and wait.
4. Confirm with a dry run on the bundled sample: `python scripts/run_headless.py`
   → it must print 5 candidates. If it does, the environment is good.

## Step 2 — Optional: AlphaMissense (calibrated missense PP3/BP4)
Only if the user wants to recover pathogenic **missense** variants (v2). It needs a
~1 GB download, so **confirm first**:
```bash
brew install htslib                         # or: conda install -c bioconda htslib
bash scripts/fetch_alphamissense.sh         # download + index (~1 GB, CC BY 4.0)
python scripts/freeze_alphamissense.py      # local, offline
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
  (If they have a GA4GH phenopacket: `python scripts/phenopacket_to_inputs.py <file>`.)

## Step 4 — Run the analysis
```bash
python scripts/run_headless.py <VCF> --hpo <HPO_FILE> --out <OUT_DIR> --timing
```
Read back the **candidate list** and per-stage timings to the user in plain terms
(gene → tier), and note the funnel (how many variants → how many candidates).

## Step 5 — Render the laudo (inline)
1. Read the generated report markdown (`<OUT_DIR>/<name>_report.md`).
2. Build a self-contained HTML report from it and publish it with the **Artifact**
   tool, using the template at
   `<repo>/.claude/skills/analyze-vcf/references/report_template.html` (from the repo
   located/cloned in Step 0) as the exact style + structure. If this skill is
   installed at the user level, that template still lives in the cloned repo.
   Fill it from the run: masthead (sample, build, generated, HPO chips), the QC
   funnel, one card per candidate (gene, coordinate, tier pill, the data-fact row,
   the `rule_path`, and the full ACMG criteria table with met / N-A / — states),
   then Methods + the disclaimer footer.
3. Give the user the Artifact link and a 2–3 line summary of the findings
   (primary vs secondary/ACMG-SF, and any ABraOM-dropped candidates).

## Guardrails (always)
- **Draft, not diagnostic.** Keep the "not for clinical use" banner in every laudo.
- **GRCh38 only.** Flag any other build.
- **Privacy.** The VCF never leaves the machine; confirm before any network step.
- Prefer running steps yourself via `Bash`; hand a command to the user only when it
  needs their data, their credentials, or a large download.
