# How vcf2report works — methods & backend

An IMRD description of the tool: what problem it solves, the local annotation and
classification backend, what it produces, and its limitations. This is the "how it works"
companion to [`BENCHMARK.md`](BENCHMARK.md) (which measures *how well* it works).

---

## Introduction

Interpreting an exome under the **ACMG/AMP** framework (Richards et al., *Genet Med* 2015) is
slow, manual, and hard to audit: a variant scientist gathers population frequencies, in-silico
predictions, ClinVar assertions, and phenotype fit, then applies ~28 criteria and a combining
rule to reach a 5-tier call. vcf2report automates the mechanical parts **deterministically and
offline**, and leaves judgment to a human — with an auditable trail for every criterion.

Three design principles shape the backend:

1. **Local-first / privacy-preserving.** The whole pipeline — parsing, QC, frequency lookup,
   in-silico scoring, ACMG classification — runs on the machine with bundled databases; no
   patient variant leaves the host unless network annotation is explicitly opted in
   (`VCF2REPORT_ALLOW_NETWORK=1`). Genomic data is sensitive; the default is offline.
2. **Never fabricate a false absence.** A variant *missing* from a frequency store is only
   treated as "absent" (which lets PM2 fire) when the store is provably complete for that
   region. An unavailable lookup returns `AF=None`, not `0.0`, so the tool cannot over-call a
   variant simply because it did not find it.
3. **Auditable, conservative, human-in-the-loop.** Every call carries the exact criterion trail
   (`rule_path`); the engine is deliberately strict; and the AI layer (Claude) *reasons* about
   free-text and narrative but **never** decides the classification.

The tool is `vcf2report` v0.1.0, GRCh38, gnomAD v4.1.

---

## Methods

### The local annotation backend

For each post-QC variant the engine attaches population frequencies, in-silico scores, clinical
assertions, and phenotype fit — every source has an offline path.

**gnomAD v4.1 — Parquet + DuckDB (the primary offline fast path).** Allele frequencies come from
a columnar **DuckDB/Parquet** store built from gnomAD v4.1 (joint exomes+genomes) by
`scripts/build_gnomad_parquet.py`. Only allele-frequency fields are extracted
(`chrom,pos,ref,alt,filter,af,af_grpmax,ac,an,nhomalt,faf95,grpmax_pop` + per-ancestry AFs),
ZSTD-compressed. At run time `gnomad_parquet.prime()` performs **one vectorised DuckDB
`LEFT JOIN`** of the entire post-QC variant set against `read_parquet(...)` — the whole exome is
annotated in a single query rather than one lookup per variant. Two safety gates matter:

- **PASS-only.** The join keeps only `filter = 'PASS'` gnomAD records, so a filtered artifact
  (AS_VQSR / AC0) can never mask a real variant's frequency (the ClinGen filtering-AF standard).
- **Region-aware absence (`mode`).** A `_meta.json` sidecar declares the store's coverage:
  `partial` (default) *never* asserts absence; `full` asserts AF=0 on any covered contig;
  `bed` asserts absence only for a variant **inside a covered exome-BED interval**. The BED is a
  kit-agnostic **MANE Select + MANE Plus Clinical** exome target (GENCODE r46, ±50 bp) built by
  `scripts/build_exome_bed.py`. Off-store variants are left unprimed → they fall through, never
  a fabricated `0.0`. This is what lets PM2 fire soundly without over-calling.

gnomAD resolves through a tiered client (`gnomad.lookup`): **Parquet → local tabix → cache →
remote tabix (public GCS bucket) → GraphQL API → bundled JSON snapshot → unknown(AF=None)**,
authoritative local stores first so a stale cache never shadows fresh data.

**AlphaMissense — PP3/BP4.** Missense pathogenicity comes from AlphaMissense (Cheng et al.,
*Science* 2023; CC BY 4.0), a tabix-indexed hg38 table; the max `am_pathogenicity` across
matching transcripts feeds **PP3** (strong ≥0.99 / moderate ≥0.90 / supporting ≥0.564) and
**BP4** (benign ≤0.34, capped at Supporting per Richards Table 5). It takes precedence over the
REVEL/CADD fallback (REVEL ≥0.70 / CADD ≥20 pathogenic; ≤0.15 / ≤10 benign). A non-missense
site returns `None`, never 0.

**ClinVar — full local store (offline) + optional live E-utilities.** Clinical assertions resolve
cache → **the full local ClinVar** (all ~4.2M GRCh38 variants, tabix-indexed for random access,
built by `scripts/build_clinvar_local.py`, ~47 MB) → live NCBI E-utilities (only if online, for a
variant newer than the snapshot) → a tiny fallback slice. The local store is complete and offline —
no network call is needed for a real ClinVar lookup — and only a positive chr+pos+ref+alt match is
accepted. Feeds the deprecated-gated **PP5** (+1 supporting) and, independently of the ACMG math,
the **ClinVar safety flag** (below).

**HPO phenotype matching — PP4 & routing.** Gene↔phenotype similarity uses the HPO `is_a` graph
with **Lin/Information-Content semantic similarity** (best-match-average) when the ontology graph
is present, falling back to exact term overlap. The average match feeds **PP4** (≥0.60); the
single strongest match routes a variant into the primary/secondary findings (≥0.50).

**ABraOM (Brazilian frequencies).** The SABE admixed-Brazilian cohort is checked alongside
gnomAD: a variant absent from gnomAD but common in Brazilians must **not** earn PM2 — a real,
local source of misclassification the tool guards against. A miss returns `None` (unknown), never
a fabricated Brazilian absence. **Note — the full ABraOM dataset is not installed yet:** the repo
ships only a 2-variant demo stub (`data/abraom/abraom_sabe.tsv`) that demonstrates the mechanism
(the OBSCN/TTN drops in the examples). Real use requires the full ABraOM SABE from IB-USP
(<http://abraom.ib.usp.br>) or `ABraOM_AF` in the VCF INFO; the mechanism is implemented, the
population coverage is not yet loaded.

**Gene constraint & pre-annotated INFO.** gnomAD LoF constraint (pLI ≥0.9 / LOEUF <0.35) drives
**PVS1**. When a VCF is already annotated (SnpEff/VEP + vcfanno), the engine reads gnomAD/ClinVar/
REVEL/CADD/AlphaMissense straight from the INFO column — a zero-lookup offline path.

### The ACMG classification engine

`classify(variant, annotation)` runs every registered criterion in a fixed, inspectable order and
combines them into a 5-tier call plus an auditable `rule_path`. Deterministic criteria are
evaluated met/not-met (PVS1, PM2, PM4, PP3, PP4, PP5, BA1, BS1, BS2, BP4, BP7); criteria that
need wet-lab or case data are *surfaced as evidence but not auto-applied* (PS1/PS3/PS4/PM1/PM5/
PP2); single-proband-inapplicable criteria (PS2/PM3/PM6, de novo / in-trans / segregation) are
marked N/A. **PVS1 strength** is modulated by the ClinGen SVI decision tree (start-loss →
Moderate; last-exon NMD-escaping → Strong; else Very Strong). Two combining models are available:
**Richards 2015 Table 5** (the conservative default) and the **ClinGen/Tavtigian 2020 points**
model (opt-in via `VCF2REPORT_ACMG_MODEL=clingen`).

**The PVS1 mechanism gate — and why population constraint alone is not it.** ClinGen SVI's first
question is whether *LoF is a known mechanism of disease for this gene*. Using pLI/LOEUF as the
proxy answers a **different question**: those metrics measure selection against **heterozygous**
LoF. In a recessive disorder the carrier is healthy, no heterozygous selection acts, and the gene
scores as tolerant — so a constraint-only gate rejects precisely the genes whose mechanism *is*
loss of function. Measured on the 100-exome cohort: 24 of 44 misses were null variants stopped at
that gate, their genes overwhelmingly recessive (Bardet-Biedl, Omenn, trichothiodystrophy...). The
same blind spot misfires on late-onset dominants — TP53's LOEUF is 0.469, i.e. "not LoF-intolerant".

So the gate is `constraint OR ClinGen Haploinsufficiency=3 OR an established autosomal-recessive
phenotype`, and the report names which route fired (`lof_mechanism_basis`):
- **ClinGen Haploinsufficiency=3** — the authoritative route: an expert panel's curated statement
  that losing one copy of the gene causes disease. This is exactly what ClinGen SVI asks PVS1 to key
  on, and it rescues the late-onset / incompletely-penetrant dominants that population constraint
  misses (TP53 has LOEUF 0.469, "not intolerant", yet is a textbook haploinsufficient tumour
  suppressor). 418 genes, refreshed by `scripts/fetch_clingen_hi.sh`.
- **population constraint** (pLI/LOEUF) — the proxy that fills gaps ClinGen has not curated.
- **established autosomal-recessive phenotype** (HPO, offline) — for recessive genes, which
  constraint is structurally blind to.

Dominant genes without HI or constraint evidence still do NOT fire: a dominant phenotype can be
gain-of-function or dominant-negative, where a null is *not* the mechanism. The constraint and HPO
routes remain gene-level proxies (a gene with both an AD gain-of-function disease and an AR LoF
disease is not yet disease-scoped — see Honest limitations); the ClinGen HI route is the curated
ground truth where it exists.

**Carrier findings.** Opening PVS1 to recessive genes makes the engine call the heterozygous null
alleles every healthy person carries (2–3 on average). The Pathogenic **tier is correct** — ACMG
classifies the variant, not the patient — but a single copy is *not* diagnostic, and phenotype
routing cannot separate it, because recessive disease genes have exactly the phenotypes a proband
presents with. Measured before the fix: a het LIPA/SKIC2 carrier was reported as the "likely
explanatory finding" while the true COPZ1 sat in "other". So a lone het P/LP in a recessive-**only**
gene is routed to its own **Carrier findings** section — kept (it has reproductive relevance),
tier untouched, but never competing with the indication. Homozygotes, a second hit in the same gene
(possible compound heterozygote), and genes with any dominant/X-linked disease are excluded, so a
real diagnosis is never hidden. This also keeps the ACMG SF v3.2 contract: the recessive SF genes
(ATP7B, MUTYH, BTD, GAA, HFE…) are reportable as actionable secondary findings **only when
biallelic** — a carrier must not be reported.

### Pipeline stages

`run_pipeline(vcf)` is an ordered, timed sequence:

```
parse → guards(build/multi-sample) → QC(DP≥10/GQ≥20/AB 0.25–0.75/PASS) →
gnomAD prime(one DuckDB join over the whole set) → annotate → filter →
AlphaMissense(candidates only) → classify(ACMG) → build report
```

The **filter funnel** keeps variants that are rare (max gnomAD/ABraOM AF ≤ 0.005) **and**
coding/splice — but **ClinVar P/LP bypasses rarity and impact** so a known pathogenic variant is
never dropped. A genome-build guard skips all coordinate-keyed annotation on a confirmed non-
GRCh38 input rather than mis-annotate; a loud warning fires if a configured gnomAD store resolves
0 of ≥50 variants (an over-calling guard).

### Local engine vs Claude — the two-phase analysis

Run from Claude (the `vcf2report` skill), the pipeline is presented as labelled phases so it is
always visible *who* performs each step:

| phase | who | what |
|---|---|---|
| Setup / liftover | 🖥️ **Local** | GRCh37→38 liftover if needed |
| Phenotype → HPO | 🤖 **Claude** | map clinician free-text to `HP:` ids (the one genuinely AI-reasoned step) |
| Frequencies | 🖥️ **Local** | detect + load the gnomAD DuckDB/Parquet store |
| Classify (ACMG) | 🖥️ **Local** | the full deterministic engine, offline, no LLM |
| Laudo | 🤖 **Claude** | synthesize the narrative report — *"do not re-classify; report the engine's calls faithfully"* |

Every number, frequency lookup, criterion evaluation, and tier call is the **deterministic local
engine**. Claude only maps free text to ontology ids and renders the final narrative. The
classification is never an LLM judgment — there are no model IDs anywhere in the engine code.

### Laudo / report generation

`build_report` assembles a `ReportModel`: findings are split into **primary** (phenotype-related,
non-benign), **secondary** (unrelated P/LP in an ACMG SF v3.2 gene — actionable, opt-in), and
**other** (incl. incidental P/LP). A deterministic `summarize()` writes the bottom-line
conclusion. An independent **ClinVar safety flag** surfaces any ≥2★ ClinVar-Pathogenic variant
the engine tiered lower — without touching the ACMG math (no circularity). `render.py` emits a
Markdown report with a **per-variant, auditable ACMG criterion table**, which the skill re-renders
into a self-contained HTML laudo. The report is labelled a *draft for expert review*.

### Tooling

Pure-Python core (only hard dependency: **Jinja2** for templating); everything else degrades
gracefully. **DuckDB** does the vectorised whole-exome frequency join and writes the Parquet in
the build scripts; **pysam/tabix (htslib)** does indexed random-access lookups for local/remote
gnomAD and AlphaMissense (imported lazily — absent → fallback); **cyvcf2** is optional fast VCF
parsing; **bcftools** is an external CLI used only by the offline build scripts.

---

## Results

Validated on **5,335 unique real pathogenic loci** (GA4GH Phenopacket Store, BSD-3) spiked into a
healthy reference exome (GIAB NA12878, public domain), with real per-locus ClinVar review status.
Reported as three honest tiers, not one number: **engine P/LP 42%**, **+ genuine ≥2★ ClinVar 52%**,
**+ phenotype 88%** (the phenotype tier is partly circular in a spike-in benchmark). Of the 1,067
genuinely ≥2★ ClinVar-known variants, **1,067/1,067 surface, 0 missed**, and the 554 the engine
tiers below P/LP are all caught by the ClinVar safety flag. Across the 2,258 PM2-driven P/LP calls,
**0** rest on a frequency above the PM2 ceiling *by the store's own reckoning* — but a whole-exome
specificity test (NA12878) reveals the store itself has false absences (see Discussion), so this
sensitivity result must be read alongside the specificity limitation. Full method, per-consequence
breakdown, and the adversarial-audit findings are in [`BENCHMARK.md`](BENCHMARK.md).

An adversarial audit of this run additionally found and **fixed** a real safety bug (a ClinVar
review-status parser that scored every real assertion 0, silently disabling the safety flag on
production inputs) and corrected an over-stated earlier headline.

---

## Discussion

**Conservative by design.** The engine's own P/LP rate (42%) is deliberately strict: single-
proband criteria (PS2/PM3/PM6/PP1/BS4) are N/A, model-adjudicated criteria are not auto-applied,
and region-aware PM2 refuses to assert absence off-panel. The gap between the strict tier and the
88% "surfaced" figure is bridged by the ClinVar safety flag and phenotype prioritization — which
is why they are reported separately.

**Honest limitations.** (1) Single-proband: no de-novo/segregation/phasing evidence. (2) The
AlphaMissense/REVEL/CADD strength thresholds are documented *seed* values pending empirical
calibration. (3) Phenotype specificity: the primary routing now uses the best-match-average, not
the single strongest match, which cut the decoy (random-phenotype) false-match from 62% to 22%
(discriminative gap 20 → 54 pts); the honest sensitivity measure remains rank against real
background (POGZ ranked #1 of 2,394 on a whole NA12878 exome). (4) **Specificity on a healthy exome:** on NA12878 the engine calls 29 P/LP
(only the planted one is real). A fixed store bug (a present-but-filtered variant read as absent →
spurious PM2) removed several (29→23); the residual are homozygous novel LoF indels in constrained
genes that are absent from *all* of gnomAD — exome-calling artifacts, addressed by upstream
variant/region-callability QC and cautious PVS1 for uncorroborated novel LoF, not by the frequency
store. A joint exomes+genomes store (`--preset joint`) still helps a separate genome-present /
exome-absent class. See
[`BENCHMARK.md`](BENCHMARK.md#specificity--phenotype-circularity--measured-on-a-real-exome).
(5) **The PVS1 mechanism gate is a proxy, and gene-level.** "The gene has an established recessive
phenotype" is strong evidence that LoF causes disease there, but it is not ClinGen gene-disease +
dosage curation, and it does not distinguish *which* of a gene's diseases is recessive — a gene with
an AD gain-of-function disease **and** an AR LoF disease currently opens the AR route on the whole
gene. The fix needs disease-scoped inheritance, which the local HPO table cannot supply (it was
deduplicated at build time and dropped `disease_id`). (6) **Late-onset dominants stay blind:** TP53
(LOEUF 0.469) is neither constraint-intolerant nor recessive, so PVS1 does not fire on a TP53 null
even though LoF is its established mechanism. (7) The gene-constraint table's **pLI column is empty**,
so `lof_intolerant` is decided by LOEUF alone. (8) Carrier routing needs zygosity from the VCF; a
sites-only VCF (no genotype) cannot be triaged this way.

**Reproducibility & privacy.** Every input is openly licensed and every store is rebuildable from
public sources (`scripts/build_*`), or downloadable as a checksummed Release. Because the default
path is fully offline, a clinical exome can be interpreted without the variant data ever leaving
the machine — the AI assists with phenotype mapping and narrative, not with moving or deciding on
patient data.
