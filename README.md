# vcf2report

**Turn a raw exome VCF into an auditable, clinical-grade ACMG variant report —
locally, offline, in seconds.**

Built with Claude — Life Sciences Hackathon.

Exome interpretation today is slow and manual: an analyst hand-filters thousands
of variants, cross-references several databases, and matches each candidate to the
patient's phenotype one by one. vcf2report turns that into minutes: point it at a
VCF and the patient's phenotype (HPO terms) and it runs the whole pipeline —
parse ▸ QC ▸ annotate ▸ prioritize ▸ classify ▸ report — and hands back a **draft
for expert review**, with the full ACMG reasoning shown for every variant.

> ⚕️ **Draft-generation aid, not a diagnostic device.** Every report is a draft to
> be verified and signed out by a qualified professional. See
> [docs/DISCLAIMERS.md](docs/DISCLAIMERS.md).

> 🤖 **Want Claude to run it for you?** Install the guided skill in **one step** and
> just say *"analyze this VCF"* — Claude sets up dependencies, runs the pipeline, and
> renders the report. See **[vcf2report.md](vcf2report.md)**.

---

## 🧬 Worked examples — real cases, real results

Six cases from the **[hpo-spiked-exomes](https://github.com/gbbarra/hpo-spiked-exomes)** benchmark
(200 real 1000 Genomes exomes, each with one known pathogenic variant planted **tell-free** + the
patient's HPO). The engine is scored **blind** — the plant carries no marker. What vcf2report returns:

| Case | Gene | Disease | Variant | vcf2report finding |
|---|---|---|---|---|
| SYN-004 | **NIPBL** | Cornelia de Lange 1 | `p.Arg1837*` (stop-gain) | 🟥 **Pathogenic** — primary diagnosis |
| SYN-073 | **BBS2** | Bardet-Biedl 2 | splice-acceptor | 🟥 **Pathogenic** — primary diagnosis |
| SYN-051 | **PIGA** | Congenital anomalies-hypotonia-seizures | missense | 🟦 Ranked #1 (phenotype-matched) — honestly held at **VUS** |
| SYN-093 | **TGFBR1** | Loeys-Dietz 1 | in-frame deletion | 🟦 Ranked #1 (phenotype-matched) — honestly held at **VUS** |
| SYN-197 | **SPINT2** | Congenital secretory diarrhea | frameshift (het) | 🟨 **Carrier** — recessive, reproductive relevance (not a diagnosis) |
| SYN-070 | **RBSN** | Congenital myelofibrosis | missense | 🟪 **VUS — flagged** for expert review (probable-pathogenic triage) |

Two confident diagnoses, two phenotype-matched candidates the engine **honestly holds at VUS**, a
recessive **carrier** it refuses to over-call, and a VUS **triaged** for review — the range, and the
restraint. Across the full 200: **177 / 200 primary recovery (88.5%)** — see
[docs/BENCHMARK.md](docs/BENCHMARK.md). Annotated example VCFs live in
[`data/example/`](data/example/); drop one into `scripts/run_headless.py` and read the laudo.

---

## What it does

Given a single-proband GRCh38 VCF + HPO terms, it produces a Markdown report with:

- a **QC & filtering funnel** (total → PASS → QC → rarity → impact → candidates);
- an **estimated sequencing-quality** panel from the VCF itself — depth & genotype
  quality at called sites, Ti/Tv, indel:SNV, % multiallelic, dbSNP novelty,
  het-allele-balance and PASS rate — a quick read on how the run sequenced;
- **primary findings** (candidates in genes overlapping the patient's phenotype);
- **secondary findings** (P/LP variants in the ACMG SF v3.2 actionable-gene list);
- a **per-variant, fully auditable ACMG trail** — all 28 ACMG/AMP criteria shown,
  each with the concrete evidence value, its source (DB + accession + date), a
  one-line reason, and whether the **engine** decided it or it is left for
  **model/expert** adjudication;
- the exact **combining-rule path** to the tier (e.g. `PVS1 + PM2 → Likely Pathogenic`).

### What makes it different

- **Auditable ACMG, not a black box.** You see *why* a variant is Pathogenic — the
  ~20 criteria determinable from a single-proband VCF are evaluated deterministically;
  the trio/segregation/phasing ones (PS2/PM3/PM6/…) are honestly marked **N/A**;
  genuine judgment calls are tagged for expert/model review.
- **Brazilian population frequencies (ABraOM).** Alongside gnomAD it checks the
  ABraOM (SABE) admixed-Brazilian cohort, so a variant *absent from gnomAD but
  common in Brazilians* is correctly dropped — the report names the spurious
  candidates a gnomAD-only pipeline would have kept.
  > ⚠️ **The full ABraOM dataset is not installed yet.** The repo ships only a
  > **2-variant demo stub** (`data/abraom/abraom_sabe.tsv`) to demonstrate the
  > mechanism (it drives the OBSCN/TTN drops in the examples). For real use, obtain
  > the full ABraOM SABE dataset from IB-USP (<http://abraom.ib.usp.br>) and populate
  > that file, or annotate the VCF with `ABraOM_AF`. A variant not in the table is
  > treated as *unknown*, never a fabricated Brazilian absence.
- **Calibrated AlphaMissense (optional).** Missense pathogenicity uses AlphaMissense
  (CC BY 4.0) at a **ClinGen-calibrated evidence strength** (PP3/BP4), validated to
  recover pathogenic missense **without** ever flipping a benign variant (below).
- **PVS1 strength by the ClinGen SVI tree.** When the VCF carries an exon rank
  (VEP `EXON` / SnpEff), null variants are graded deterministically — Very Strong for
  NMD-triggering, **downgraded to Strong** for last-exon (NMD-escaping) nonsense/
  frameshift and **Moderate** for start-loss (Abou Tayoun 2018). Un-annotated VCFs
  stay Very Strong, so nothing is silently over-called.
- **Ontology-aware phenotype matching (HPO).** PP4 and the primary-vs-secondary
  routing use a Lin/Information-Content similarity over the HPO `is_a` graph, so a
  patient term matches a *related* gene term (parent/child) weighted by specificity —
  and adding more (explained) phenotypes no longer dilutes the score, the failure
  mode of plain term overlap on phenotype-rich cases. Falls back to exact overlap
  when the graph isn't built.
- **Offline gnomAD, no 150 GB download.** A local **DuckDB/Parquet** store of gnomAD
  v4.1 frequencies (29.6M variants) resolves a whole exome's frequencies in *one*
  vectorised join — a real exome classifies **fully offline in ~7 s**. Get it two ways,
  same result, into the auto-detected local `data/gnomad/gnomad_parquet/`:
  `scripts/fetch_gnomad_parquet.sh` (download a checksummed copy) **or**
  `scripts/build_gnomad_parquet.py` (rebuild from the public bucket — prove it on one
  chromosome in minutes). A reduced-tabix path (`build_gnomad_local.py`) and GRCh37
  liftover (`liftover_to_grch38.py`) also exist. A safety model makes a false-absence
  impossible (partial stores never assert absence; a configured-but-missing store warns
  loudly instead of silently over-calling). See
  [docs/LOCAL_ANNOTATION.md](docs/LOCAL_ANNOTATION.md) · gnomAD ODbL-1.0:
  [data/gnomad/NOTICE.md](data/gnomad/NOTICE.md).
- **Private by default.** Runs fully offline on bundled data; the VCF never leaves
  the machine, and outbound lookups are opt-in (`VCF2REPORT_ALLOW_NETWORK=1`).
- **Auditable & tested.** 246 automated tests; the classification logic is validated
  against real ClinVar ground truth (below) and was hardened by multi-agent
  adversarial review.

---

## Observed performance (hackathon)

**Concordance panel** — the engine run against **200 real ClinVar variants**
(100 pathogenic / 100 benign), with ClinVar *withheld* from the engine so the
comparison is non-circular (details: [docs/CONCORDANCE.md](docs/CONCORDANCE.md)):

| Metric | Result | Meaning |
|---|---|---|
| **Gross discordances (P↔B flips)** | **0 / 200** | never calls a benign variant pathogenic, or vice-versa |
| **Pathogenic precision** | **100%** | when it says Pathogenic-ish, it is |
| **Benign precision** | **100%** | when it says Benign-ish, it is |
| **Concordance when decisive** | **100%** | of the calls it commits to, all match ClinVar |
| **Pathogenic sensitivity** | **60%** (LoF 84%) | conservative — defers the rest to VUS |
| v1 → v2 (AlphaMissense) | 37% → **60%** | +23 missense recovered, still 0 gross |

The engine is deliberately **conservative** (it defers to VUS unless it has strong
deterministic evidence) — the value is that it is **never dangerously wrong** and
is **right whenever it commits**. Speed: the bundled sample runs end-to-end in
**well under a second** (~0.4 s on a laptop). All numbers are reproducible:
`python3 scripts/run_concordance.py`.

**See a real rendered laudo** (self-contained HTML — open in a browser):
[SYN-001](docs/example_reports/SYN-001.report.html) ·
[SYN-002](docs/example_reports/SYN-002.report.html). More
[example reports](docs/example_reports/) (Markdown) are bundled too.

---

## Install

You need **Python ≥ 3.10**. Pick the tier you want — the demo needs only the first.

```bash
git clone https://github.com/gbbarra/vcf2report.git
cd vcf2report
python3 -m pip install -e .            # core engine (only dependency: jinja2)
```

That is enough to run the offline demo below — **all databases for it are already
bundled** in `data/` (a tiny ClinVar fallback slice, gnomAD snapshot, ABraOM, HPO,
gene constraint, in-silico; ~2 MB total, no download).

### Data stores — REQUIRED for real-exome analysis

A **real exome** needs the three full annotation Parquet stores; the `/vcf2report` store gate
**blocks an analysis until they are present + intact**. Get all three with **one command** after
installing:

```bash
bash scripts/setup_stores.sh
```

| Store | Size | Source |
|---|---|---|
| **gnomAD v4.1** | ~1.3 GB | pre-built GitHub release — `fetch_gnomad_parquet.sh` (checksummed, frozen) |
| **ClinVar GRCh38** | ~60 MB | pre-built GitHub release — `fetch_clinvar_parquet.sh` (checksummed, rebuilt **weekly**) |
| **AlphaMissense hg38** | ~1 GB | fetched from DeepMind + built locally — **CC BY-NC-SA 4.0, not redistributed** |

Prereqs: `gh` (GitHub CLI), `zstd`, and `duckdb` (`pip install duckdb`); AlphaMissense also needs
htslib/tabix (`brew install htslib` / `conda install -c bioconda htslib`). Check each store's
availability, version, build date, and integrity anytime with `python3 scripts/check_stores.py`.

Optional extras:

```bash
python3 -m pip install -e ".[fast]"    # + cyvcf2 (htslib-backed VCF reader, faster)
python3 -m pip install -e ".[mcp]"     # + MCP SDK (to drive it from Claude Desktop)
python3 -m pip install -e ".[dev]"     # + pytest/hypothesis (to run the test suite)
```

### Optional databases (for richer / real-exome analysis)

| Database | Needed for | How |
|---|---|---|
| **Live gnomAD / ClinVar** | up-to-the-minute frequencies | opt in with `VCF2REPORT_ALLOW_NETWORK=1` (sends only variant coordinates, never the VCF) |
| **SnpEff + reference** | annotating a *raw* (un-annotated) VCF | `conda install -c bioconda snpeff bcftools vcfanno htslib` — see [docs/LOCAL_ANNOTATION.md](docs/LOCAL_ANNOTATION.md) |

If your VCF is already annotated (SnpEff/VEP `ANN`/`CSQ`, gnomAD/AlphaMissense in
`INFO`), vcf2report reads those directly and needs none of the above.

---

## Run it

### 1. The bundled demo (offline, zero setup)

```bash
python3 scripts/run_headless.py
```

```
Report written to data/out/sample_exome_report.md
  candidates classified: 5
  - SCN1A: Pathogenic
  - PAX6: Likely Pathogenic
  - RB1: Likely Pathogenic
  - KCNQ2: Uncertain Significance (VUS)
  - CACNA1A: Uncertain Significance (VUS)
```

Open `data/out/sample_exome_report.md` to read the full auditable report.

### 2. Your own VCF + phenotype

```bash
python3 scripts/run_headless.py path/to/exome.vcf --hpo hpo_terms.txt --out out/
python3 scripts/run_headless.py path/to/exome.vcf --hpo hpo_terms.txt --stdout   # print to screen
```

- **VCF**: single-proband, **GRCh38**. A multi-sample VCF analyses the first column
  (or pass the proband's sample name); a non-GRCh38 build is flagged and its
  coordinate lookups are skipped.
- **`hpo_terms.txt`**: one HPO id per line (free text after the id is ignored):

  ```
  HP:0001250   Seizure
  HP:0002133   Status epilepticus
  ```

  (Have a GA4GH Phenopacket instead? `python3 scripts/phenopacket_to_inputs.py`.)

Flags: `--hpo FILE` · `--out DIR` · `--stdout` · `--sample-id ID` · `--timing`.

### 3. From Claude — two ways

**a) Claude Code (recommended — one-step install, then plain language).** Install the
guided skill once; it works in any session and bootstraps everything itself:

```bash
mkdir -p ~/.claude/skills/vcf2report && curl -fsSL \
  https://raw.githubusercontent.com/gbbarra/vcf2report/main/.claude/skills/vcf2report/SKILL.md \
  -o ~/.claude/skills/vcf2report/SKILL.md
```

Then say *"analyze this VCF: /path/to/exome.vcf"* — Claude clones/installs (if
needed), runs the pipeline, and **renders the laudo inline**. See **[vcf2report.md](vcf2report.md)**.

**b) Claude Desktop (natural-language chat via MCP).** `pip install -e ".[mcp]"`, add the
`mcpServers` block from [`claude_desktop_config.example.json`](claude_desktop_config.example.json),
restart, and ask in plain language. Full guide: [docs/SETUP.md](docs/SETUP.md).

---

## Architecture

Three front doors, one engine:

```
Claude Code   ── /vcf2report skill (guided harness) ──┐
Claude Desktop ── Agent Skills ▸ MCP server ───────────┤
Terminal      ── scripts/run_headless.py ──────────────┤
                                                       ▼
                          vcf2report package
                          parse ▸ qc ▸ annotate ▸ filter ▸ ACMG ▸ report
                          ├─ live APIs (opt-in): gnomAD, ClinVar, HPO
                          └─ local: ClinVar slice, gnomAD cache, ABraOM,
                             HPO, gene constraint, AlphaMissense
```

All logic lives in the importable, headless-testable `vcf2report` Python package;
the Claude Code skill and the MCP server are thin adapters over the same pipeline.
Details: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

**Standards:** ACMG/AMP (Richards et al., Genet Med 2015) + ClinGen SVI refinements
· ACMG SF v3.2 (Miller et al., 2023) · ClinGen PP3/BP4 calibration (2024) · HGVS ·
GA4GH Phenopackets.

## Documentation

| Doc | What |
|---|---|
| [vcf2report.md](vcf2report.md) | **Run it from Claude** — one-step skill install + the guided harness |
| [SETUP.md](docs/SETUP.md) | Full install + Claude Desktop / MCP integration |
| [DEMO_SCRIPT.md](docs/DEMO_SCRIPT.md) | Guided demo walkthrough |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | How the pieces fit together |
| [CONCORDANCE.md](docs/CONCORDANCE.md) | The ClinVar validation panel + AlphaMissense calibration |
| [LOCAL_ANNOTATION.md](docs/LOCAL_ANNOTATION.md) | Annotating a raw VCF with SnpEff |
| [ANNOTATION.md](docs/ANNOTATION.md) | Local vs live annotation trade-offs |
| [PHENOPACKET.md](docs/PHENOPACKET.md) | GA4GH Phenopacket input |
| [SYNTHETIC_CASES.md](docs/SYNTHETIC_CASES.md) | How the bundled synthetic exomes were built |
| [DISCLAIMERS.md](docs/DISCLAIMERS.md) | Scope & limitations |
| [example_reports/](docs/example_reports/) | Rendered (HTML) + Markdown example laudos |

## Test it

```bash
python3 -m pip install -e ".[dev]"
python3 -m pytest -q        # 246 tests
```

## License

MIT (see [LICENSE](LICENSE)). Bundled/optional datasets keep their own licenses
(AlphaMissense CC BY 4.0; ClinVar/gnomAD/HPO public; ABraOM per its terms).
