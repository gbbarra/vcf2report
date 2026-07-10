# Run vcf2report *with Claude* (guided harness)

vcf2report is designed to be driven **from Claude** — you don't memorise commands.

## Install into Claude — one step

Install the skill once (copies a single file into Claude's user skills folder — it
then works in **any** Claude Code session, in any folder):

```bash
mkdir -p ~/.claude/skills/analyze-vcf && curl -fsSL \
  https://raw.githubusercontent.com/gbbarra/vcf2report/main/.claude/skills/analyze-vcf/SKILL.md \
  -o ~/.claude/skills/analyze-vcf/SKILL.md
```

Restart Claude Code. That's it — you don't need to clone the repo or install
anything else yet; the skill does the rest on first use.

> Already have the repo open in Claude Code? You can skip the install — the same
> skill ships in this repo's `.claude/skills/` and is picked up automatically.

## Use it — just ask

In any Claude Code session, type `/analyze-vcf` or just say, in plain language:

> *"Analyze this VCF for a patient with seizures: `/path/to/exome.vcf`"*
> *"Set up vcf2report and run the demo."*
> *"Gera um laudo desse exoma: `~/caso.vcf.gz`, HPO HP:0001250, HP:0011097."*

Claude picks up the **`analyze-vcf`** skill (`.claude/skills/analyze-vcf/SKILL.md`)
and acts as a guided harness — it sets up the environment, runs the pipeline, and
renders an auditable ACMG report (laudo) inline.

## What it looks like

> **You** — *"analyze this VCF for a patient with seizures: `~/SYN-001.vcf.gz`
> (HPO HP:0001250, HP:0002133, HP:0011097)"*
>
> **Claude** — checks Python + the package (installs if missing), writes the HPO
> file, then runs:
> ```
> python3 scripts/run_headless.py ~/SYN-001.vcf.gz --hpo /tmp/hpo.txt --out out/ --timing
> ```
> *"Filtered **24,801 variants → 2 candidates** in 2.8 s:*
> - ***SCN1A*** *(frameshift) → **Pathogenic** — matches the seizure phenotype (PVS1 + PM2 + PP4 + PP5)*
> - ***RB1*** *(frameshift) → **Pathogenic** — actionable secondary finding (ACMG SF v3.2)"*
>
> …and renders the full auditable laudo inline → 🧬 **a clinical report** with the
> QC funnel, findings by category, and every variant's ACMG criterion trail
> ([see a rendered example](docs/example_reports/SYN-001.report.html)).

You never touched a command reference — you described the case, Claude did the rest.

## What Claude does for you (the harness flow)

1. **Environment & install** — checks Python ≥ 3.10 and the package; if missing it
   hands you the exact `git clone` + `pip install -e .` (only dependency: jinja2;
   the demo databases are bundled — no download). It confirms with a dry run.
2. **Optional databases** — offers AlphaMissense (calibrated missense, ~1 GB) and
   live gnomAD/ClinVar, both **opt-in** and confirmed before anything downloads or
   touches the network. Your VCF never leaves the machine.
3. **Inputs** — asks for your **GRCh38** VCF and the patient's **HPO terms** (or maps
   your free-text phenotype to HPO ids; a GA4GH phenopacket works too).
4. **Run** — executes `scripts/run_headless.py` and reads the candidates + funnel
   back to you in plain terms.
5. **Laudo** — renders a clean, auditable HTML report inline (per candidate: gene,
   variant, ACMG tier, and the full criterion trail with met / N-A / — states),
   using `.claude/skills/analyze-vcf/references/report_template.html`.

Anything that needs your machine, your data, or a big download, Claude gives you the
exact command to run in your terminal and waits for the output — everything else it
runs itself.

## Prefer the terminal directly?

```bash
git clone https://github.com/gbbarra/vcf2report.git && cd vcf2report
python3 -m pip install -e .
python3 scripts/run_headless.py                                   # bundled demo
python3 scripts/run_headless.py my.vcf --hpo hpo.txt --out out/   # your case
```

See the full [README](README.md) for install tiers, the observed performance, and
the Claude Desktop / MCP integration ([docs/SETUP.md](docs/SETUP.md)).

> ⚕️ vcf2report produces a **draft for expert review — not a diagnostic device.**
