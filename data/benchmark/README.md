# Benchmark baseline — `hpo-spiked-exomes`

`hpo-spiked-exomes-baseline.tsv` is the **per-case result** behind the `179/200` recorded in
[`docs/BENCHMARK.md`](../../docs/BENCHMARK.md). 17 KB, one row per case, checked in so a future
run can be diffed against it instead of against a number in prose:

```bash
python3 scripts/run_benchmark.py \
    --annotated <bench>/realistic_annotated --bench <bench> --jobs 4 \
    --out after.tsv --compare data/benchmark/hpo-spiked-exomes-baseline.tsv
```

`--compare` prints net primary recovery before→after and names every case that moved bucket or
tier (★ marks missense). A change of headline count tells you *that* something moved; this file is
what tells you **which**.

## What this run was

| | |
|---|---|
| engine | `main` @ `0ebee07` (identical result on `7dcc39d`, pre-#24 — 0 cases differ) |
| cohort | [gbbarra/hpo-spiked-exomes](https://github.com/gbbarra/hpo-spiked-exomes) release `data-v1` (2026-07-22), SHA256-verified |
| stores | gnomAD v4.1 69,898,057 · AlphaMissense hg38 71,034,269 · ClinVar 4,195,020 (built 2026-07-14 — **18 d old**, past the 14 d policy) |
| result | primary 179 · carrier 1 · probable-VUS 5 · other 8 · absent 7 · **errors 0** |

`withheld_tier` is the tier the planted variant gets with its **own** ClinVar assertion nulled
(`--withhold-clinvar`) — the novel-variant scenario. It is identical to `tier` in every row, which
is the measurement retiring the claim that PP5 carries these plants; see `docs/BENCHMARK.md`.

## Reading `outcome` honestly

`outcome` records **where the planted gene landed**, not whether it was called pathogenic — 70 of
the 179 `primary` rows carry `tier = Uncertain Significance (VUS)`. Filter on both columns to get
classification accuracy (109/200). And 45 of the 200 plants have an answer-key ClinVar record that
is not itself P/LP, one of them `Benign/Likely_benign`, so 200 is not a reachable ceiling.
