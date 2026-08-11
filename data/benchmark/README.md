# Benchmark baseline — `hpo-spiked-exomes`

`hpo-spiked-exomes-baseline.tsv` is the **per-case result** behind the `180/200` recorded in
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
| engine | branch head after the coordinate-matched scorer; the engine tree is unchanged from `ac1c6bf` |
| criteria off | **PS1 / PM1 / PM5 could not fire** — the genome-wide ClinVar residue index was never built here (NCBI unreachable). Three of 28 criteria disabled is a condition of this number. |
| cohort | [gbbarra/hpo-spiked-exomes](https://github.com/gbbarra/hpo-spiked-exomes) release `data-v1` (2026-07-22), SHA256-verified |
| stores | gnomAD v4.1 69,898,057 · AlphaMissense hg38 71,034,269 · ClinVar 4,195,020 (built 2026-07-14 — **18 d old**, past the 14 d policy) |
| result | primary 180 · carrier 1 · probable-VUS 5 · other 8 · absent 6 · **errors 0** |

`withheld_tier` is the tier the planted variant gets with its **own** ClinVar assertion nulled
(`--withhold-clinvar`) — the novel-variant scenario. It equals `tier` in **193 of 195** scored rows,
which is the measurement retiring the claim that PP5 carries these plants; see `docs/BENCHMARK.md`.
The two exceptions (SYN-172 *CHD7*, SYN-197 *SPINT2*) drop Pathogenic → Likely Pathogenic, i.e. PP5
was the one supporting line that completed a Pathogenic rule path — it still never decided P/LP
versus not.

## Reading `outcome` honestly

`outcome` records **where the planted gene landed**, not whether it was called pathogenic — 77 of
the 180 `primary` rows carry `tier = Uncertain Significance (VUS)`. Filter on both columns to get
classification accuracy (103/200). And 45 of the 200 plants have an answer-key ClinVar record that
is not itself P/LP, one of them `Benign/Likely_benign`, so 200 is not a reachable ceiling.

## `outcome` and `tier` answer different questions

They are deliberately not the same key, because the cohort makes them different:

- **`outcome`** is **gene-level** — did the planted gene reach the primary bucket? That is the
  published metric. A clinician shown the right gene has the diagnosis to confirm, whichever allele
  carried it there.
- **`tier` / `consequence`** describe the **planted variant**, matched on `(chrom, pos, ref, alt)`
  from `manifest/planted_variants.tsv`. **40 of the 200 cases plant two alleles** (compound het —
  the manifest's `allele` column says `primary`/`second`), so "the classification for this gene" is
  a question the cohort gave two answers to on purpose.
- **`tier_source`** says which: `planted-locus` (193 rows), `gene-fallback` (2 — the coordinate was
  not classified, so the gene-level pick stands in), or empty (5 — nothing matched at all).

Before this, the bucket came from one variant and the tier from another; SYN-001 *TM2D3* was
recorded `primary·Pathogenic` when what it has is a Pathogenic variant in `carrier` and a separate
VUS in `primary` — a pair present in no single row.

That coordinate match also exposed a scoring bug: a row can only carry `outcome=absent` *and* a
tier if the variant was classified under a **different gene symbol**. Two did. One is a real
annotation limitation (SYN-042: the key plants in *RNU4-2*, the annotator attributes the locus to
the adjacent paralogue *RNU4-1*). The other, SYN-182, was `C10ORF71` in the key versus HGNC's
`C10orf71` from the annotator — an exact string comparison, and a recovered frameshift sitting in
the primary bucket scored as a miss. Matching is case-insensitive now; that is the +1 from 179.
