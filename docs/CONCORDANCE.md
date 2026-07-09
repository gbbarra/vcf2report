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

## What the panel measures

With ClinVar withheld, the engine's power lives on a few axes:

- **LoF pathogenicity** — a null variant in a LoF-intolerant gene earns PVS1, and
  rare frequency earns PM2 → Likely Pathogenic / Pathogenic.
- **Common-benign** — a variant common in gnomAD earns BA1 / BS1 → Benign /
  Likely Benign.
- **Missense pathogenicity (v2, AlphaMissense)** — see below.

### v1 (gnomAD-only) and the missense gap

Without an in-silico score, a **pathogenic missense** variant earns only PM2 and
lands in **VUS** — the engine defers missense pathogenicity to evidence it does
not have. On the frozen panel this shows as a high *LoF* pathogenic sensitivity
but a low *overall* one, reported honestly rather than hidden.

### v2: AlphaMissense with ClinGen-calibrated PP3/BP4

v2 adds **AlphaMissense** (CC BY 4.0; `annotate/alphamissense.py`, local tabix on
the hg38 file) as the missense in-silico axis, applied at a **variable ACMG
strength** per the ClinGen 2024 recalibration: a strong pathogenic prediction can
reach **PP3_Strong**, which — combined with PM2 (rare) — lifts a missense from VUS
to **Likely Pathogenic**. The score→strength cutoffs live in `config.py`
(`AM_PP3_STRONG` etc.) as documented **seed** values: they are meant to be tuned
against this panel (raise them until gross discordances stay at zero) and verified
against the ClinGen table before clinical use. AlphaMissense benign evidence (BP4)
is capped at Supporting, since Richards Table 5 has no benign-moderate bucket.

The panel is the guardrail for this calibration: if the thresholds were too loose,
a rare ClinVar-**benign** variant with a high AlphaMissense score would be elevated
to Pathogenic — a **gross discordance** the test asserts stays at zero.

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

### v2: freeze AlphaMissense (offline, local file)

The missense in-silico axis needs the AlphaMissense hg38 file downloaded once
(~1 GB), then a purely local freeze:

```bash
bash scripts/fetch_alphamissense.sh          # download + tabix-index (~1 GB, CC BY 4.0)
python scripts/freeze_alphamissense.py        # local tabix -> alphamissense_frozen.json (no network)
```

`data/concordance/alphamissense_frozen.json` is small and committed; the ~1 GB
source file is gitignored. If it is absent, the panel simply runs v1
(frequency-only) behaviour.

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
| `src/vcf2report/annotate/alphamissense.py` | AlphaMissense client (local tabix) for the v2 missense axis |
| `scripts/build_concordance_panel.py` | One-time freeze (network): ClinVar harvest + gnomAD freeze |
| `scripts/fetch_alphamissense.sh` / `freeze_alphamissense.py` | Download + offline freeze of AlphaMissense scores |
| `scripts/run_concordance.py` | Offline renderer |
| `tests/test_concordance_panel.py`, `tests/test_alphamissense_v2.py` | Harness math, engine path, calibrated PP3/BP4, and the frozen-panel safety guard |
| `data/concordance/` | `ground_truth.tsv` + `gnomad_frozen.json` + `alphamissense_frozen.json` (generated once, then committed) |
