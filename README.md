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

---

## What it does

Given a single-proband GRCh38 VCF + HPO terms, it produces a Markdown report with:

- a **QC & filtering funnel** (total → PASS → QC → rarity → impact → candidates);
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
- **Calibrated AlphaMissense (optional).** Missense pathogenicity uses AlphaMissense
  (CC BY 4.0) at a **ClinGen-calibrated evidence strength** (PP3/BP4), validated to
  recover pathogenic missense **without** ever flipping a benign variant (below).
- **Private by default.** Runs fully offline on bundled data; the VCF never leaves
  the machine, and outbound lookups are opt-in (`VCF2REPORT_ALLOW_NETWORK=1`).
- **Auditable & tested.** 138 automated tests; the classification logic is validated
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
`python scripts/run_concordance.py`.

An example report lives in [docs/example_reports/](docs/example_reports/).

---

## Install

You need **Python ≥ 3.10**. Pick the tier you want — the demo needs only the first.

```bash
git clone https://github.com/gbbarra/vcf2report.git
cd vcf2report
python -m pip install -e .            # core engine (only dependency: jinja2)
```

That is enough to run the offline demo below — **all databases for it are already
bundled** in `data/` (ClinVar slice, gnomAD snapshot, ABraOM, HPO, gene constraint,
in-silico; ~2 MB total, no download).

Optional extras:

```bash
python -m pip install -e ".[fast]"    # + cyvcf2 (htslib-backed VCF reader, faster)
python -m pip install -e ".[mcp]"     # + MCP SDK (to drive it from Claude Desktop)
python -m pip install -e ".[dev]"     # + pytest/hypothesis (to run the test suite)
```

### Optional databases (for richer / real-exome analysis)

| Database | Needed for | How |
|---|---|---|
| **AlphaMissense hg38** (~1 GB, CC BY 4.0) | calibrated missense PP3/BP4 (v2) | `bash scripts/fetch_alphamissense.sh` then `python scripts/freeze_alphamissense.py` (needs htslib: `brew install htslib` / `conda install -c bioconda htslib`) |
| **Live gnomAD / ClinVar** | up-to-the-minute frequencies | opt in with `VCF2REPORT_ALLOW_NETWORK=1` (sends only variant coordinates, never the VCF) |
| **SnpEff + reference** | annotating a *raw* (un-annotated) VCF | `conda install -c bioconda snpeff bcftools vcfanno htslib` — see [docs/LOCAL_ANNOTATION.md](docs/LOCAL_ANNOTATION.md) |

If your VCF is already annotated (SnpEff/VEP `ANN`/`CSQ`, gnomAD/AlphaMissense in
`INFO`), vcf2report reads those directly and needs none of the above.

---

## Run it

### 1. The bundled demo (offline, zero setup)

```bash
python scripts/run_headless.py
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
python scripts/run_headless.py path/to/exome.vcf --hpo hpo_terms.txt --out out/
python scripts/run_headless.py path/to/exome.vcf --hpo hpo_terms.txt --stdout   # print to screen
```

- **VCF**: single-proband, **GRCh38**. A multi-sample VCF analyses the first column
  (or pass the proband's sample name); a non-GRCh38 build is flagged and its
  coordinate lookups are skipped.
- **`hpo_terms.txt`**: one HPO id per line (free text after the id is ignored):

  ```
  HP:0001250   Seizure
  HP:0002133   Status epilepticus
  ```

  (Have a GA4GH Phenopacket instead? `python scripts/phenopacket_to_inputs.py`.)

Flags: `--hpo FILE` · `--out DIR` · `--stdout` · `--sample-id ID` · `--timing`.

### 3. From Claude Desktop (natural language)

`pip install -e ".[mcp]"`, add the `mcpServers` block from
[`claude_desktop_config.example.json`](claude_desktop_config.example.json) to your
Claude Desktop config, restart, and say: *"Analyze this VCF for a patient with
seizures: /path/to/exome.vcf"*. Full guide: [docs/SETUP.md](docs/SETUP.md).

---

## Architecture

```
Bench scientist (Claude Desktop)  ── natural language + VCF path ─┐
                                                                  ▼
Agent Skills (the clinical SOP) ──►  MCP server (thin adapter) ──► vcf2report package
                                                                   parse ▸ qc ▸ annotate
                                                                   ▸ filter ▸ ACMG ▸ report
                                                                   ├─ live APIs: gnomAD, ClinVar, HPO
                                                                   └─ local: ClinVar slice, gnomAD cache,
                                                                      ABraOM, HPO, constraint, AlphaMissense
```

All logic lives in the importable, headless-testable `vcf2report` Python package;
the MCP server is a thin adapter. Details: [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md).

**Standards:** ACMG/AMP (Richards et al., Genet Med 2015) + ClinGen SVI refinements
· ACMG SF v3.2 (Miller et al., 2023) · ClinGen PP3/BP4 calibration (2024) · HGVS ·
GA4GH Phenopackets.

## Documentation

| Doc | What |
|---|---|
| [SETUP.md](docs/SETUP.md) | Full install + Claude Desktop / MCP integration |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | How the pieces fit together |
| [CONCORDANCE.md](docs/CONCORDANCE.md) | The ClinVar validation panel + AlphaMissense calibration |
| [LOCAL_ANNOTATION.md](docs/LOCAL_ANNOTATION.md) | Annotating a raw VCF with SnpEff |
| [ANNOTATION.md](docs/ANNOTATION.md) | Local vs live annotation trade-offs |
| [PHENOPACKET.md](docs/PHENOPACKET.md) | GA4GH Phenopacket input |
| [SYNTHETIC_CASES.md](docs/SYNTHETIC_CASES.md) | How the bundled synthetic exomes were built |
| [DISCLAIMERS.md](docs/DISCLAIMERS.md) | Scope & limitations |

## Test it

```bash
python -m pip install -e ".[dev]"
python -m pytest -q        # 138 tests
```

## License

MIT (see [LICENSE](LICENSE)). Bundled/optional datasets keep their own licenses
(AlphaMissense CC BY 4.0; ClinVar/gnomAD/HPO public; ABraOM per its terms).
