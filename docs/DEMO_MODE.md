# Demo mode — showing the pipeline without the stores

The guided flow's Stage-1 gate is hard: it refuses to analyse anything unless the gnomAD,
AlphaMissense and ClinVar Parquet stores are all present and intact. That is deliberate — a laudo
produced with gnomAD absent has **PM2 / BA1 / BS1 / BS2** blind, and a blind rarity criterion
over-calls rather than under-calls.

But the gate cannot tell a patient exome from `data/example/SYN-073.BBS2.annotated.vcf.gz` — a file
committed to this repository purely so the pipeline can be demonstrated. So on any machine without
~1 GB of stores (a fresh container, a cloud session driven from a phone) the demo cases were
unrunnable through the guided flow. That is the one place a demo is most useful.

**Demo mode** is the single, narrow exemption.

## Using it

```bash
# the guided flow: pass demo: true in the analysis Workflow's args
# the gate, directly:
python3 scripts/check_stores.py --gate --demo data/example/SYN-073.BBS2.annotated.vcf.gz

# the CLI:
python3 -m vcf2report.cli data/example/SYN-073.BBS2.annotated.vcf.gz \
    --hpo data/example/SYN-073.hpo.txt --demo --out out/demo
```

`VCF2REPORT_DEMO=1` is the environment equivalent, for surfaces that cannot pass a flag.

The available fixtures and their planted diagnoses are listed in
[SYNTHETIC_CASES.md](SYNTHETIC_CASES.md). `SYN-073` (BBS2, Bardet-Biedl) reaches a Pathogenic call
and also exercises carrier triage; `SYN-070` (RBSN) lands at VUS and exercises the
probable-pathogenic triage queue.

## The two properties that make it safe

Both are enforced in `src/vcf2report/demo.py` and covered by `tests/test_demo_mode.py`, not left to
convention.

### 1 · It can only ever apply to a committed fixture

`--demo` is **not a store override**. Pointed at anything outside `data/example/` it is *refused*,
not honoured — the gate exits non-zero and the CLI exits 2:

```
⛔ demo mode REFUSED — '/data/patient.vcf.gz' is not one of this repository's committed example
   VCFs. Demo mode is not a store override: it exists only to demonstrate the pipeline on its own
   fixtures. Build the Parquet stores to analyse this file.
```

Paths are `resolve()`d before the check, so a relative path, a symlink, or `data/example/../..`
cannot smuggle a file in from outside.

So the flag can never be the reason a real exome was analysed with blind stores. Someone who
believed they were in demo mode is told plainly that they are not, rather than silently getting an
ordinary run.

### 2 · The stamp is unconditional, and does not depend on the flag

Any run whose input is a `data/example/` VCF is stamped — **whether or not `--demo` was passed**,
whether or not the stores happen to be present, and on every surface (CLI, headless script, MCP).
The provenance is derived in `pipeline.run_pipeline` from the VCF that was actually parsed and the
stores that were actually found, then carried on `ReportModel.provenance`.

This asymmetry is the point. A forgotten flag has to fail toward *stamped*; the dangerous direction
is an unstamped fixture laudo passing for a patient result. A genuine analysis mislabelled "demo"
merely under-claims.

The stamp reaches the reader in five places, each independently tested:

| Where | What it says |
|---|---|
| Conclusion, **first bullet** | `DEMONSTRATION RUN — NOT A PATIENT RESULT` + which criteria are unverified |
| Masthead banner, `templates/report.md.j2` | the same, above the findings |
| Masthead banner, the built-in renderer | the same (this is what MCP and `--stdout` emit) |
| Methods → `data_mode` + two `criteria_from_*` keys | `DEMONSTRATION — … stores absent (…) — nothing looked up live` |
| `<sample>_results.json` → `provenance` | the machine-readable record, with the two groups kept apart |
| The Artifact, `{{DEMO_BANNER}}` | a full-width warning, deliberately louder than the DRAFT chip |

## What the stamp names — and the distinction that matters

A vague "demo data" note would leave a reader guessing which rows of the trail to distrust. But
"frozen demo slices" is *also* wrong, for most of it. With the stores absent, a demo run's evidence
comes from **two different places**, and they fail differently:

### Baked into the fixture's own INFO fields

`scripts/make_annotated_example.py` was run on a machine where the stores **were** present. It read
the **complete** databases and wrote the values into the example VCF's INFO — the exact fields the
engine reads back through its `from_vcf` path:

```
gnomad_AF  gnomad_AC  gnomad_AN  gnomad_nhomalt  gnomad_faf95
CLNSIG  CLNREVSTAT  CLNVI  CLNDN
am_pathogenicity
```

So these are **real measurements from the full databases** — e.g. the BBS2 diagnostic variant carries
`gnomad_AF=8.99282e-07` against `gnomad_AN=1461862`, the genuine gnomAD v4.1 allele number. What they
are *not* is re-checkable: they are frozen at generation time, cover only this fixture's variants, and
carry no version or build date you can verify on a store-less machine.

| Store absent | Criteria resting on baked INFO |
|---|---|
| gnomAD | PM2, BA1, BS1, BS2 |
| AlphaMissense | PP3, BP4 |
| ClinVar (assertion) | PP5, BP6 |

### From the frozen ClinVar residue slice

The residue index behind **PS1 / PM1 / PM5** is *not* carried in the VCF — it is a `(gene, aa_pos,
alt_aa)` lookup. The repo ships `data/clinvar/clinvar_residue_frozen.tsv.gz`, ~3.6k rows covering
**15 demo genes** (APC, BBS2, FBN1, KCNQ2, NIPBL, PIGA, RB1, SCN1A, SCN2A, SLC2A1, SPINT2, STK11,
STXBP1, TGFBR1, WT1). For any other gene these three report *"index unavailable — not assessed"* —
honestly not-assessed, never a negative.

This is the one part that genuinely is a demo slice.

---

With all three stores absent that is **11 of the 28 criteria**, split 8 baked / 3 sliced. The
provenance block carries them apart (`criteria_from_fixture_info`, `criteria_from_frozen_slice`),
Methods gets one key each, and `tests/test_demo_mode.py` asserts the two lists never overlap — a code
in both would mean the split is wrong.

Stage 5 (analysis triage) reports all 11 as **UNVERIFIED**, never as "available".
