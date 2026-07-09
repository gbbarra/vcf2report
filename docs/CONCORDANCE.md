# Concordance panel — ClinVar ground truth vs the engine

The concordance panel is the **quantitative** counterpart to vcf2report's
per-variant audit trail. It takes a set of real, expert-classified ClinVar
variants and asks: *when the engine sees each one, does it reach the same broad
call — and does it ever get one dangerously wrong?*

## Why it is not circular

The ACMG engine consumes ClinVar directly through **PP5** (a supporting line that
fires when ClinVar already reports Pathogenic). Scoring the engine against ClinVar
*while feeding it ClinVar* would inflate agreement for free. So the panel runs
every variant with ClinVar **withheld** from the annotation — the engine must
recover the call from population frequency (gnomAD/ABraOM), LoF mechanics and
in-silico evidence alone. ClinVar is used only as the answer key.

## What v1 measures (and what it deliberately does not)

v1 is **gnomAD-only** (no bundled in-silico scores). With ClinVar withheld, the
engine's *deterministic* power lives on two axes:

- **LoF pathogenicity** — a null variant in a LoF-intolerant gene earns PVS1, and
  rare frequency earns PM2 → Likely Pathogenic / Pathogenic.
- **Common-benign** — a variant common in gnomAD earns BA1 / BS1 → Benign /
  Likely Benign.

A **pathogenic missense** variant, with ClinVar withheld and no in-silico score,
earns only PM2 and lands in **VUS**. That is by design: the engine defers missense
pathogenicity to in-silico / hotspot evidence (`adjudicated_by="model"`) it does
not have here. The panel reports this honestly rather than hiding it — expect a
high *LoF* pathogenic sensitivity and a lower *overall* pathogenic sensitivity.

## The safety invariant

The one number that must always hold: **zero gross discordances**. A gross
discordance is a PATH↔BEN flip — the engine calling a ClinVar-benign variant
Pathogenic, or a ClinVar-pathogenic variant Benign. A VUS on either side is a
conservative non-call; a flip is a real error. `tests/test_concordance_panel.py`
asserts this is zero on the frozen panel.

## Building the panel (once, with network)

The build harvests real ClinVar variants (coordinates and labels **always from
ClinVar — never fabricated**) for a curated gene panel, then freezes each
variant's gnomAD grpmax frequency via remote tabix. It writes two files under
`data/concordance/`:

```bash
VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_concordance_panel.py
# ~100 gnomAD lookups, ~1-2 min. Idempotent/resumable — re-run to top up.
#   data/concordance/ground_truth.tsv    (the answer key)
#   data/concordance/gnomad_frozen.json  (frozen frequencies)
```

Both files are committed. After this one-time step, everything below is **fully
offline, forever**.

## Running the panel (offline)

```bash
python scripts/run_concordance.py                      # Markdown matrix + metrics
python scripts/run_concordance.py --json               # machine-readable
python scripts/run_concordance.py --out data/out/concordance.md
python scripts/run_concordance.py --include-clinvar    # circular; for comparison only
```

The offline guard test runs with the suite:

```bash
pytest tests/test_concordance_panel.py                 # SKIPs until the panel is frozen
```

## Files

| File | Role |
|---|---|
| `src/vcf2report/concordance.py` | Harness: load panel, classify (ClinVar withheld), confusion matrix + metrics, Markdown |
| `scripts/build_concordance_panel.py` | One-time freeze (network): ClinVar harvest + gnomAD freeze |
| `scripts/run_concordance.py` | Offline renderer |
| `tests/test_concordance_panel.py` | Harness math, engine path, builder parsing, and the frozen-panel safety guard |
| `data/concordance/` | `ground_truth.tsv` + `gnomad_frozen.json` (generated once, then committed) |
