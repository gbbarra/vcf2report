# How the validation cohort was built — an honest methods record

This document records **exactly** how the SYN validation exomes were constructed, so the
benchmark numbers (`docs/BENCHMARK.md`) can be read for what they are. The cohort is
**synthetic**: real, healthy exome backgrounds with a **spiked-in** pathogenic variant. It is
**not** real patient data, and every figure derived from it should be read with the construction
below in mind.

> TL;DR — the background exomes are **real 1000 Genomes samples downloaded (streamed) from a public
> S3 bucket**; they are **not** derived from the gnomAD/ClinVar Parquet stores (those are the
> engine's frequency/clinical references, a different thing). One pathogenic variant, taken from a
> real published patient case, is then **inserted** into each background.

## 1. Background exome — real, downloaded, never repeated

Each case starts from a **different** real human exome: a **1000 Genomes DRAGEN v4.4.7** sample.
The caller's `hard-filtered.vcf.gz` (~415 MB per sample) is **streamed over HTTPS from the public
S3 bucket** `1000genomes-dragen-v4-4-7` — no `aws` CLI, no credentials, no reference FASTA needed.
The stream is resumable with stall-detection (abort < 1 MB/s for 25 s, retry with `-C -` on a fresh
connection); a sample that stays too slow is skipped and retried on a later pass, so the whole run
survives S3's throttled windows (`scripts/make_syn_cohort.sh` + `run_cohort_v2_loop.sh`).

Every case uses a **distinct** sample — 100 distinct backgrounds in v1, another 100 distinct ones in
the v2 expansion (SYN-101..200), none reused. This is deliberate: a planted variant is tested against
**100 different real genetic backgrounds**, not one shared background, so background-specific quirks
cannot inflate the result.

**This is a download, not a Parquet derivation.** The gnomAD v4.1 and ClinVar Parquet stores the
engine reads are population-frequency / clinical-assertion references — they are never turned into
the cohort VCFs. The cohort VCFs are genuine per-sample caller output.

## 2. Subset to the engine's region — one streaming pass

`bcftools norm -m -any` (split multiallelics) piped into `bcftools view -T data/gnomad/exome_hg38.bed`
subsets each ~4–5 M-variant genome VCF to the **vendor-neutral MANE / GENCODE v46 exome BED** — the
same MANE Select + MANE Plus Clinical ±50 bp region the engine covers. Result: ~100 k variants per
case, a realistic exome scale. No kit-specific (e.g. IDT/Twist) BED is used anywhere.

## 3. The spike — one exact pathogenic variant from a real patient case

The causative variant is drawn from a **real GA4GH Phenopacket-Store case** (release 0.1.27, ~10 k
published patient cases). For each chosen gene we take:

- the **exact coordinate** `chrom:pos:ref:alt` the phenopacket records;
- that case's **HPO terms** (the phenotype the report is given);
- the gene, molecular consequence, and disease name.

`scripts/spike_variant.py` then **inserts that exact record** into the background exome and looks the
coordinate up in the **ClinVar GRCh38 VCF** (auto-downloaded from NCBI, public domain):

- **coordinate present in ClinVar** → the spike carries ClinVar's real `CLNSIG` / `CLNREVSTAT` /
  `CLNDN`, so PP5 fires on the genuine assertion and the real disease name reaches the report;
- **coordinate absent from this ClinVar release** → the variant is still planted, but with a
  **synthetic `CLNSIG=Pathogenic`** label. This is a construction artifact and is disclosed here:
  those cases do not carry an independent ClinVar assertion.

The record is flagged `SPIKED=1` (v2 second alleles `SPIKED2=1`) in the VCF INFO — the plant is never
hidden. The sample column is de-identified to the `SYN-NNN` id.

**What this means for the benchmark:** the planted variant is a *constructed* record, so the cohort
tests the engine's **classification + prioritisation logic against a known answer**, not variant
*calling*. To keep the phenotype metric honest, an anti-circularity control runs every planted case
against a **decoy** (random, mismatched) phenotype (`docs/BENCHMARK.md`), and the "engine-only" rate
excludes the ClinVar read-back.

**Is the visible `SPIKED` marker a "tell" the engine exploits? No — verified, not assumed.** The
engine reads none of the plant's markers (`SPIKED`, `GENE`, `CSQ`, `CLNSIG`, …): it reclassifies from
coordinate + genotype + the SnpEff consequence, joining gnomAD / ClinVar / AlphaMissense / HPO from
its **own local stores**. Strip *every* marker from a planted record — leaving only `CHROM/POS/REF/ALT`,
the genotype, and the SnpEff `ANN` any annotated VCF carries — and the classification is
**byte-identical**: same tier, same criteria, same evidence values. (The ClinVar the engine reports is
its store's, not the VCF's — proven because stripping the VCF's `CLNSIG` leaves it unchanged; the
`gnomAD AF` was never in the VCF at all.) The `SPIKED=1` flag is a **build-time** convenience for
locating the plant; correctness is scored by **coordinate** from the truth TSV, never by the flag. For
a cohort with no visible marker, `scripts/spike_variant.py --realistic` plants a **tell-free** record
that borrows a real background call's INFO/FORMAT — indistinguishable from a genuine call, statistically
as well as visually. It changes no number, precisely because the markers were already ignored.

## 4. Faithful genotype (v2) — restoring the real zygosity

v1 planted every variant as a lone **heterozygote**. That silently misrepresents recessive cases: a
single het in an autosomal-recessive gene is a healthy **carrier**, not a diagnosis. v2 fixes this by
matching each planted variant back to its phenopacket case and restoring the patient's **real
genotype** (`scripts/build_v2_biallelic.py` + `data/synthetic_cohort/v2_faithful_plan_*.json`):

- **compound-heterozygous** — both real alleles planted (two het records);
- **homozygous** — the real `allelicState`, planted `1/1`;
- **single-allele** — the source genuinely recorded one allele.

The v2 distribution and the diagnostic-sensitivity trajectory it produced (38 → 58 → 75 → 91/100) are
in `docs/BENCHMARK.md`.

## 5. Variant selection — stratified, not random

The v1 100 were chosen to span molecular consequences (35 missense · 25 stop-gained · 25 frameshift ·
10 in-frame indel · 5 start-loss = 65 SNV / 35 indel) — a deliberate spread, not a random draw. The
v2 expansion **oversamples the VUS-producing classes** (67/100 missense + in-frame vs 45 in v1,
`scripts/select_cohort_v2.py`) to build up N on the cases the engine conservatively defers on.

## 6. Reproduce it

```bash
# v1 (SYN-001..100):
bash scripts/make_syn_cohort.sh            # streams S3, subsets, spikes → data/synthetic_cohort/
# or fetch the prebuilt corpus + checksum:
bash scripts/fetch_syn_cohort.sh

# v2 faithful/expanded:
bash scripts/run_cohort_v2_loop.sh         # SYN-101..200 (detached, resumable)
bash scripts/build_v2_biallelic.py --plan data/synthetic_cohort/v2_faithful_plan_101_200.json --out <dir>
bash scripts/fetch_syn_cohort_v2.sh        # or fetch the prebuilt v2 corpus
```

## 7. Honesty caveats (read before quoting any number)

- **Synthetic, not real patients.** Real backgrounds, but the causative variant is inserted. De-
  identified; **not for clinical use**.
- **Tests classification, not calling.** A spiked record is perfectly called by construction, so the
  cohort says nothing about upstream variant-calling sensitivity.
- **The default plant is visually identifiable** (it carries `SPIKED`/`GENE`/`CSQ`/`CLN*` INFO and a
  minimal FORMAT), but the engine **provably ignores** all of it — strip the markers and the call is
  byte-identical (§3). `spike_variant.py --realistic` removes them entirely for a plant indistinguishable
  from a real background call, without changing any result.
- **Some spikes carry a synthetic ClinVar label** (coordinate not in the ClinVar release) — disclosed
  in §3; the engine's non-ClinVar ("engine-only") metric is reported separately for that reason.
- **Public sources only:** 1000G DRAGEN (public S3), GA4GH Phenopacket Store (open), ClinVar (public
  domain), HPO, gnomAD/AlphaMissense (their own licenses). No private or patient-identifiable data.
