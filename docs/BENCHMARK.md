# Benchmark — does vcf2report flag real pathogenic variants?

A validation harness that answers the only question that matters for a variant-interpretation
tool: **when a real pathogenic variant is present, does the tool bring it to the clinician's
attention?** It uses real, published cases — not hand-invented variants.

## Data (all openly licensed)

> **How the cohort was built — read [`COHORT_CONSTRUCTION.md`](COHORT_CONSTRUCTION.md) first.** The
> backgrounds are real 1000G exomes **downloaded** from public S3 (not derived from the Parquet
> stores); one pathogenic variant from a real patient case is **spiked in**. The construction — and
> its honesty caveats (synthetic spike, exact-coord ClinVar vs synthetic label, tests classification
> not calling) — is documented there so these numbers are read for what they are.

| piece | source | license |
|---|---|---|
| variants + phenotypes | **GA4GH Phenopacket Store** (10,377 curated cases from 959 publications; GRCh38 coordinates, HGVS, HPO terms, ACMG P/LP) | BSD-3 |
| healthy background exome | **GIAB NA12878 / HG001** (NIST) | public domain |
| ClinVar status | ClinVar GRCh38 VCF (NCBI) | public domain |
| phenotype ontology | HPO `genes_to_phenotype` | open |

Each synthetic case = the NA12878 exome + **one** real pathogenic variant spiked in (with a
constructed SnpEff `ANN`) + that case's HPO terms. Half the cases are fed a ClinVar
Pathogenic assertion ("with ClinVar"), half are not ("without" → the engine must classify
from mechanism: PVS1/PM2/PP3/PP4). All variants are absent from the local gnomAD store (they
are rare disease alleles) — verified: on the 12-case set, 11/12 are absent from **all** of
gnomAD and 0/12 are genomes-only.

## The metric: *surfaced to the clinician*, not the raw ACMG tier

vcf2report classifies **independently** and conservatively. Region-aware PM2 fires only for a
variant genuinely absent inside a covered exome-BED interval (never a false absence off-panel),
so the raw tier stays deliberately strict. The clinically-meaningful metric is whether the
variant is **brought to attention** by *any* of:

- the engine's ACMG tier is Pathogenic / Likely Pathogenic, **or**
- it is flagged as a ≥2-star **ClinVar Pathogenic** (surfaced in the conclusion), **or**
- it is a **phenotype-matched** primary candidate (HPO overlap ≥ threshold, not benign).

**Honesty caveat on the three signals.** They measure different things and must be reported
separately. Engine P/LP is the pure ACMG call — fully engine-independent. The ClinVar flag is
credited **only for genuinely ≥2-star assertions** (star count read from the real review
status; see the fix below). The phenotype signal is the softest: in this spike-in benchmark the
HPO terms come from the same publication that names the causative gene, so an HPO overlap is
near-guaranteed by construction — phenotype-surfacing here demonstrates prioritization, not
classification. The headline therefore reports all three tiers, not a single conflated number.

## Adversarial review

Two adversarial passes: a 23-agent review of the 12-case run, then a 12-agent audit of the full
5335-case run (8 tail auditors + code auditors for harness integrity, over-calls, and the ClinVar
surface). The tail audit found **0 real misses**. The code audit found — and this pass **fixed** —
a real production bug plus two harness-fidelity issues that had inflated an earlier headline:

- **FIXED — `clinvar_stars` scored every real assertion 0 (clinical-safety bug).** The star
  function matched only underscore tokens, but both production paths deliver a *space*-delimited
  review status (the VCF-INFO path normalizes with `.replace("_"," ")`; live E-utilities returns
  spaces natively). So on real input the ≥2-star ClinVar safety flag **never fired**. Fixed to
  normalize first (as PP5 already did); regression-tested against both delimiter forms. Impact on
  the benchmark: of the 1067 genuinely ≥2-star ClinVar-P variants, the 554 the engine tiers below
  P/LP are now flagged **554/554** (previously 0).
- **CORRECTED — harness no longer fabricates a uniform 2-star assertion.** The earlier harness
  injected "2-star Pathogenic" for every ClinVar-matched locus; the real review status is 0-star
  for 1020, 1-star for 978, ≥2-star for only 1067 of 3065. The benchmark now feeds the **real**
  per-locus significance + review status, so the ClinVar credit is honest.
- **DISCLOSED — phenotype signal is partly circular** (see the honesty caveat above): reported as
  a separate tier, not folded into an engine-independent number.
- **FIXED — ClinVar surface (clinical-safety).** A ≥2-star ClinVar Pathogenic variant the
  engine tiered VUS was being reported as "no Pathogenic finding". The conclusion now flags
  it explicitly ("⚠️ Classified Pathogenic in ClinVar … DO NOT dismiss") without touching the
  ACMG math (no circularity). This alone lifted *surfaced* from **0/12 → 9/12**.
- **APPLIED — region-aware PM2 (`mode=bed`).** The exomes+MANE store was `mode=partial`, so
  PM2 could never fire on genuine absence — the tool's primary use case. It now asserts
  absence (→ PM2) only for a variant **inside a covered panel-BED interval** (where the store
  is complete for gnomAD exomes); off-panel it stays unprimed. On the benchmark this lifted
  engine **P/LP 0 → 6/12** and **surfaced 9 → 11/12** — the LoF cases now reach LP via
  PVS1 + PM2. The concordance panel is unchanged (it uses the frozen gnomAD, not the store):
  still 0 gross discordances / 100% precision / 60% sensitivity. Residual risk (a variant in
  gnomAD genomes but not exomes) was 0/12 empirically; a joint (exomes+genomes) MANE build
  closes it fully.
- **MEASURED, kept opt-in — ClinGen points model.** With region-aware PM2 in place,
  `VCF2REPORT_ACMG_MODEL=clingen` gives the *same* benchmark result (6/12 P/LP, 11/12
  surfaced) but **regresses the concordance panel to 37%** (a no-phenotype artifact: strict
  ClinGen points call a rare missense-without-phenotype VUS). So Richards stays the default;
  clingen remains available via the env var for labs that prefer the SVI points scheme.

## Results

**12-case set (mixed with/without ClinVar):** found 12/12 · engine P/LP **6/12** · **surfaced
11/12** (after the ClinVar surface + region-aware PM2). The 1 residual is a non-coding-RNA
gene (out of the protein-coding ACMG scope). Progression as the fixes landed:
0/12 surfaced → 9/12 (ClinVar surface) → 11/12 (region-aware PM2).

**Full set — 5335 unique loci** (from 9595 hg38 P/LP-or-causative records with ≥3 HPO terms;
42 SV/CNV/symbolic/>50 bp records were explicitly dropped and counted, never reclassified;
4260 duplicate loci collapsed). 677 genes; every consequence class incl. start-loss and
in-frame indels. Because a spiked variant's ACMG call is independent of the background, the
harness batch-primes gnomAD + AlphaMissense once and classifies each variant in isolation —
**5335 cases in ~18 s**. ClinVar significance + review status are fed **as they really are**,
per-locus, in the space-delimited production form.

Reported as three honest tiers, not one conflated number:

| signal | surfaced | note |
|---|---|---|
| **engine P/LP** (pure ACMG) | **2266 / 5335 (42%)** | fully engine-independent |
| **+ genuine ≥2-star ClinVar** | **2820 / 5335 (52%)** | ClinVar credited only for real ≥2-star assertions |
| **+ phenotype match** | **4725 / 5335 (88%)** | phenotype tier is partly circular here (see caveat) |

**The clinical-safety property, honest and demonstrated.** Of the 1067 genuinely ≥2-star
ClinVar-Pathogenic variants, **1067/1067 surface — 0 missed**. The 554 of them the engine
tiers below P/LP are caught by the ClinVar safety flag **554/554** (before the `clinvar_stars`
fix: 0/554). This is a real property on real review-status data, no longer a harness artifact.

By consequence (surfaced incl. phenotype): stop-gain 760/821 (92%) · frameshift 1195/1306
(91%) · missense 2605/3003 (86%) · in-frame 129/162 (79%) · start-loss 36/43 (83%) — no class
fails systematically.

**Over-calls (false positives).** Across all 2258 PM2-driven P/LP calls, **none** rests on a
population frequency above the PM2 ceiling (store check, 0 violations); an earlier independent
gnomAD-API spot-check corroborated genuine absence (bulk API re-verification is rate-limited).
The 610 not-surfaced cases are engine-VUS without a ≥2-star ClinVar assertion or phenotype
overlap — the conservative tail, not misses.

## Specificity & phenotype circularity — measured on a real exome

The isolated-classification harness measures sensitivity (does a real pathogenic variant surface)
but not specificity. Two negative controls quantify what it can't:

**Whole real exome (GIAB NA12878 + one planted variant).** We spiked a known pathogenic
POGZ variant (White-Sutton syndrome) into the real GIAB NA12878 exome (28,566 variants) with the
case's phenotype, and ran the full pipeline (~10 s). The engine **ranked POGZ #1 of 2,394
candidates** (Pathogenic, phenotype 1.00) — the honest, non-circular sensitivity signal is this
**rank against real background**, not a match rate. But it also called **29 P/LP on this healthy
individual** (only POGZ is real). We traced every one:

- **A real store bug — FIXED.** The PASS-gate turned a *present-but-filtered* variant into a false
  absence: e.g. SON was called Pathogenic on a 99.98%-frequency AS_VQSR record the store read as
  absent → spurious PM2. The join now matches regardless of filter and only asserts absence when
  *no* record exists (PASS records still definitive; non-PASS → AF unavailable, not absent). This
  removed those over-calls (**29 → 23 P/LP**) with no benchmark regression (genuine ≥2★ ClinVar-known
  held at 1067/1067). See `annotate/gnomad_parquet.py`.
- **Residual — calling artifacts, now flagged.** The rest are **homozygous** LoF indels in constrained
  genes (CYFIP2, LTBP4, PRKDC…) that a wide remote-tabix scan confirms are **absent from all of gnomAD**
  (76k genomes) — the classic signature of exome-calling artifacts in hard regions. A homozygote requires
  the allele to exist in the population, so AC=0 + hom is implausible for a real allele. **Fix (implemented):**
  such variants keep their ACMG tier in the ranked table but are routed out of the confident
  "likely-explanatory" findings into a *"verify the genotype (orthogonal/Sanger) first"* caution — on the
  POGZ demo the likely-explanatory list drops from 8 genes to 3. Heterozygous variants (incl. genuine novel
  dominant LoF like the planted POGZ) are untouched. The few remaining P/LP are *heterozygous* novel LoF,
  which are genotype-indistinguishable from a true positive and need region-callability QC to separate.
- **A joint (exomes+genomes) store** (`build_gnomad_parquet.py --preset joint`, the script default)
  is still recommended: it removes the *other* false-absence class — a variant genuinely present in
  gnomAD genomes but not exomes. The shipped store is exomes-only today (`data/gnomad/NOTICE.md`).

**Phenotype decoy control — measured, then fixed.** For all 5,335 cases we re-scored each with a
*random, unrelated* phenotype. Originally the primary-findings routing used the *single strongest*
patient↔gene term match (`hpo_best_match`), which a decoy cleared **62%** of the time (vs 83% true
— only ~20 pts specific), because one broad term matches almost any gene. **Fix (implemented):**
route on the best-match-**average** (`hpo_match_score`, requiring the phenotype *as a whole* to fit
the gene) at 0.6 — the same bar as PP4. This drops the decoy false-match to **22%** (77% true), so
the discriminative gap nearly triples, **20 → 54 pts**. The engine's own P/LP rate and the genuine
≥2-star ClinVar surface remain the fully-specific signals, reported separately.

## Whole-exome scale test — 100 planted variants on NA12878

The isolated harness above measures classification; this measures **real-world surfacing at scale**.
100 different real pathogenic variants (100 unique genes, GA4GH Phenopacket Store) were each planted
one-at-a-time into the real GIAB NA12878 exome (~28,600 variants) with that case's phenotype, and run
through the full offline pipeline (~6 s each). The honest question: is the planted variant brought to
the clinician's attention — by ACMG tier (Likely Pathogenic / Pathogenic), or as a phenotype-matched
primary candidate ranked for curation against the ~2,400 real background candidates?

**73/100 surfaced** — 46 as P/LP by tier, 8 flagged as a ≥2-star ClinVar-Pathogenic the engine tiered
lower (from the full local ClinVar), and 19 as phenotype-matched primary candidates; 27 not surfaced.
When it surfaces, the planted variant sits at **median rank 3** of the ~2,400 real background candidates
(top of the list in 27/100) — the specialist reads it near the top, not buried.

| consequence | surfaced | why |
|---|---|---|
| stop-gain | 22/25 (88%) | LoF → PVS1 + PM2 reaches P/LP |
| in-frame indel | 8/10 (80%) | PM4 + PM2 / phenotype (see the fix below) |
| start-loss | 4/5 (80%) | |
| missense | 24/35 (69%) | AlphaMissense PP3 + phenotype |
| frameshift | 15/25 (60%) | LoF, but PVS1 needs the gene flagged LoF-intolerant in the constraint table |

The misses are disclosed, not hidden:
- **Specificity, not a regression (77 → 73).** An earlier run surfaced 77/100. Tightening the
  phenotype routing to require the best-match **average** (not the single strongest term) to clear
  threshold removed **4 non-specific phenotype matches** — the intended trade: on a decoy (random,
  unrelated) phenotype the false-match rate falls **62% → 22%**. The 4 dropped were promiscuous
  single-term hits, not genuine gene↔phenotype fits; the ACMG-tier and ClinVar signals (54/100) are
  unchanged.
- **In-frame indels — fixed.** An earlier run scored 0/10 here; the cause was a consequence-term
  mismatch, not an ACMG gap. The impact filter and PM4 recognised only the VEP terms
  `inframe_insertion`/`inframe_deletion`, so an in-frame indel spelled `inframe_indel` (or SnpEff's
  `disruptive_/conservative_inframe_*`) was dropped before classification. Both now match any term
  containing "inframe" → **8/10**.
- **LoF that doesn't reach P/LP**: a nonsense/frameshift in a gene not flagged LoF-intolerant gets PM2
  only → VUS; phenotype misses track gene↔HPO annotation coverage. Both are data-coverage gaps, not
  classification errors, and a specialist curates the ranked list regardless.

This whole-exome, real-background figure (73%, rank against real candidates) is the honest, non-circular
counterpart to the isolated harness's phenotype tier — the metric that matters is that the true variant
is *brought forward* near the top, with a clear, measured breakdown of where it is not.

## 100 distinct backgrounds, real annotation — the SYN cohort

The NA12878 test above spikes into a single background with a constructed `ANN`. This one removes both
crutches: 100 **distinct** 1000-Genomes DRAGEN exomes (one planted variant each, exact ClinVar
coordinate, distinct gene + phenotype), **functionally annotated end-to-end with SnpEff on the MANE
transcript set** (`GRCh38.mane.1.5.refseq`, the same transcripts the gnomAD store is sliced by) — 99.99%
of ~100k background variants per exome annotated, so the engine can (and must) classify the whole
background, not just the spike. This is what makes specificity measurable at all: an unannotated
background cannot be called P/LP, so any "no over-call" claim on one is an artifact, not a result.

**Read the two questions separately — conflating them inflates the number.**

*Diagnostic sensitivity* — is the planted variant presented as the likely **diagnosis** (primary section)?
**91/100** on the faithful cohort. Where the planted variant lands: **91 primary · 2 incidental P/LP · 7 miss**
(plus recessive carriers routed to the carrier section on the by-design single-allele cases).

**Fidelity matters — and the number earned its way up honestly.** Each planted variant was matched back to
its GA4GH Phenopacket-Store case to recover the patient's *real* genotype: **22 compound-heterozygous** (both
true alleles), **35 homozygous** (the real `allelicState`), **43 single-allele** (the source genuinely
recorded one). An earlier cohort planted every variant as a lone heterozygote, which silently flattened
every recessive patient into a *carrier* — the engine correctly declined to call those a diagnosis, but the
cohort couldn't test recessive diagnostic recovery at all. The trajectory, all measured:

| cohort / fix | diagnostic (primary) |
|---|---:|
| lone-het spike, stale HPO table | 38/100 |
| + current HPO release (gene→phenotype + inheritance) | 58/100 |
| + **faithful biallelic genotypes** (phenopacket-matched) | 75/100 |
| + **homozygous-diagnosis routing fix** | **91/100** |

By genotype: **compound-het 22/22**, **homozygous 32/35**. All **9 non-primary are honest limitations, not
misses the engine should have made** (each verified): 2 non-coding snRNA (RNU5B-1/RNU4-2 — outside the
engine's protein-consequence scope), 2 genes the current HPO release *dropped* (U2AF2/MAPK8IP3 — still
surfaced as P/LP, just not phenotype-elevated), 4 missense **VUS** (RBSN/MOCS3/ELFN1/SETD2 — insufficient
variant-level evidence; homozygosity does not rescue a VUS), and 1 phenotype scoring just under the 0.6
routing threshold (TRAF7 at 0.572 — the measured efficient frontier, below). This is the anti-circular
counterpart to the isolated harness: no ClinVar read-back, no single shared background, and the phenotype
earns its keep — a decoy (random) phenotype routes the gene to primary only **18%** of the time vs the true
phenotype's **37%**.

*The homozygous-diagnosis fix.* Making the AR cases faithfully homozygous exposed a real interaction bug: the
`is_hom_absent_artifact` QC guard — hom + gnomAD-absent = an implausible genotype / calling-artifact
signature, demoted out of primary — was created for a **healthy** exome, but that is *also* the textbook
signature of a recessive **diagnosis** in an affected proband. It buried 16 of 35 homozygous diagnoses in
"other". The fix keeps the guard but lets a **phenotype-matched P/LP** hom-absent variant into primary (with
the confirm-the-genotype caveat still attached); a hom-absent variant that does *not* match the phenotype
stays demoted, so calling artifacts still cannot flood a healthy proband's report.

*What else moved the number.* Two data/logic fixes, each measured, not assumed:
- **PVS1 recognises recessive-LoF disease.** Gating PVS1 on population constraint (pLI/LOEUF) is blind to
  recessive genes — the carrier is healthy, so the gene never looks constrained (see HOW_IT_WORKS). Adding
  an established-AR-phenotype route lets the engine classify recessive LoF as the pathogenic variant it is
  (then route lone hets to carrier). This is what produced the 36 correct carrier calls.
- **PVS1 keys on ClinGen Haploinsufficiency=3** — the curated "LoF causes disease" statement PVS1 actually
  asks for, added as a third mechanism route. It rescues the late-onset / incompletely-penetrant dominants
  constraint misses (TP53: LOEUF 0.469, yet HI=3). On the concordance panel this lifted pathogenic
  sensitivity 63% → 67% (LoF-only 90.9% → **100%**) with **zero** pathogenic↔benign flips and unchanged
  benign precision. Measured directly on healthy exomes too: P/LP-per-exome median stayed **4** (max 10),
  and a het LoF landing in an HI gene — the new incidental class — occurred once in 15 exomes (TGFBR1),
  so opening 418 dominant genes did not flood the background.
- **The gene→phenotype/inheritance table was stale.** It was a frozen pyhpo snapshot missing 36 of the
  cohort's disease genes entirely — dropping both their phenotype match and their inheritance. Rebuilding
  from the current HPO release covered all 36 (each with inheritance): diagnostic sensitivity 38 → 58,
  misses 27 → 4, decoy specificity essentially flat (17% → 18%) while true-phenotype routing rose 29% → 37%.

*Phenotype aggregation is already near-optimal — measured, not assumed.* The primary routing scores a gene
by the **best-match average** over the patient's HPO terms, at a 0.6 threshold. Sweeping alternatives on
the cohort (top-K averages, max) against the decoy control refutes the tempting "top-K" idea: it lifts
sensitivity but lifts the decoy false-match *more* (neuro/developmental terms are promiscuous, so a random
phenotype also has a few strong matches). The all-term average is what separates a true fit (matches many
of the patient's terms) from a decoy (matches a few) — it sits on the efficient frontier, and if anything
0.65 buys more specificity (decoy 32% → 21% for one lost case).

**Over-call (specificity), now real:** median **4** P/LP per ~100k-variant healthy exome (max 10) — of
which ~60% are heterozygous recessive **carrier** alleles, which every healthy person carries a few of
(normal biology, correctly tiered and routed to carrier, not diagnosis). The prior "median 1, no
flooding" was a pure artifact of the 0%-annotated background and is not comparable.

*Scale check — a second, fully independent 100.* The cohort was doubled to 200: SYN-101..200 are **100
more distinct 1000G backgrounds** (none reused), **100 new genes**, variants **oversampled toward the
VUS-producing consequences** (67/100 missense + in-frame, vs 45 before) to stress the cases the engine
defers on. Diagnostic sensitivity on this independent set is **89/100** — statistically the same as the
first 100's 91, so the number is a property of the engine, not a fluke of one cohort (**~180/200**
overall). Every non-primary case is again an honest limitation, not a miss the engine should have made.
(Building the new 100 also earned its keep as a stress test: the new backgrounds hit a ClinVar record
whose disease-name field carried a stray newline, which split the spiked VCF record in **14** cases — a
cohort-construction bug the first 100 never triggered; found, fixed at the root, and re-hardened.)

## Probable-pathogenic VUS — triaged, not reclassified

The conservative engine correctly leaves a variant at **Uncertain Significance** when its evidence is all
*Supporting* (a strong AlphaMissense score, a phenotype match, and even a ClinVar Pathogenic assertion do
not, individually, reach Likely-Pathogenic — and PP5 is capped for anti-circularity). But a VUS carrying
several such signals is not the same as one carrying none. `report/vus_triage.py` ranks the VUS by that
suggestive evidence, **itemised** (so the operator sees *why*, not a black-box score), and surfaces the
top ones for expert + model exploration — **the ACMG tier is never changed**. Two gates keep it specific
to the indication: the gene must be **phenotype-relevant**, and it must carry **molecular support** beyond
that — which drops the flag rate from 8.9 to **0.9 per exome**. On the oversampled new 100, the planted
missense variant that the engine held at VUS is correctly prioritised in **4** cases (54/56 missense
surfaced as diagnosis-or-triaged). This is the tool's stated stance made operational: the deterministic
engine does the ACMG call; the genuinely uncertain variants are handed forward — with their evidence — for
the judgement the engine deliberately withholds.

## Honest limitations

- Single-proband: PS2/PM3/PM6/PP1/BS4 (de novo, in-trans, segregation) are N/A.
- In-frame indels rely on PM4 (Moderate) — no PVS1 — so without a phenotype match or ClinVar they can
  stay VUS; a repeat/hotspot-aware PM4 refinement would help further (now 8/10 after the term fix above).
- Specificity on a healthy exome is limited (above) — disclosed and quantified, not hidden.
- The exomes+MANE store keeps PM2 conservative (see above) — a deliberate no-false-absence
  choice, quantified here rather than hidden.
- **The SYN cohort plants lone heterozygotes**, so it cannot test diagnostic recovery for recessive
  genes (a lone het there is a carrier by definition). 36/100 cases are affected; a cohort v2 must plant
  AR cases biallelically. The engine's carrier routing is correct — the cohort is what is incomplete.
- **The PVS1 recessive-LoF route is a gene-level proxy** (see HOW_IT_WORKS): "the gene has an established
  recessive phenotype" is not disease-scoped, so a gene with both an AD gain-of-function disease and an AR
  LoF disease opens the route on the whole gene. Disease-scoped inheritance needs an HPO table rebuild that
  keeps `disease_id` (the local one dropped it). TP53-style late-onset dominants stay a constraint blind
  spot, and the constraint table's pLI column is empty (LOEUF only).
- The NA12878 harness above still uses a constructed `ANN`; the SYN cohort is the SnpEff-MANE-annotated
  measurement that supersedes it (and confirms the exon rank it omitted correctly downgrades start-loss
  PVS1 — measured, not assumed).

## Reproduce

The harness (survey phenopackets → curate a balanced set → spike into NA12878 → run + score)
is deterministic given the three public inputs above; see `scripts/` and this doc.
