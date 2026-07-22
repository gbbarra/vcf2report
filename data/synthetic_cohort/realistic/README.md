# Synthetic exome benchmark — HPO-linked pathogenic variants spiked into real backgrounds

A ground-truth corpus for **validating variant annotators, variant callers, and clinical-report
generators**. Each sample is a **real, healthy 1000 Genomes exome** with **one known pathogenic
variant** (and, where the source case is biallelic, its real second allele) **planted at an exact
coordinate**, carrying that case's **HPO phenotype**. The answer key is known; in the *realistic*
build the plant is **not marked** in the VCF, so a tool can be scored **blind**.

> **Synthetic, not real patients.** Real public backgrounds + inserted variants. De-identified.
> **Not for clinical use.**

---

## What's in a release

| Path | What it is |
|---|---|
| `realistic/SYN-NNN.vcf.gz` | **Raw, tell-free** VCF — the plant carries a real DRAGEN call's INFO/FORMAT and **no marker** (indistinguishable from a genuine call). Un-annotated. |
| `realistic_annotated/SYN-NNN.annotated.vcf.gz` | The **same** VCF **annotated with SnpEff** (`GRCh38.mane.1.5.refseq`) — adds `ANN/LOF/NMD`, still tell-free, DRAGEN INFO/FORMAT preserved. |
| `marked/SYN-NNN.vcf.gz` | The **original marker-bearing** VCF (the plant carries `SPIKED=1;GENE=;CSQ=;CLN*` INFO). Kept for transparency / easy inspection. |
| `planted_variants.tsv` | **The truth manifest** — every planted allele: `chrom:pos:ref:alt`, gene, zygosity, primary/second, consequence, ClinVar significance/review/id. |
| `sidecars/SYN-NNN.planted.tsv` | Per-sample truth (same columns). |
| `sidecars/SYN-NNN.hpo.txt` | Per-sample HPO terms (one `HP:` per line). |
| `cohort*.tsv` | Per-sample config (sample id, gene, coord, consequence, disease, HPO). |

The raw is **byte-derivable** from the annotated (`bcftools annotate -x INFO/ANN,INFO/LOF,INFO/NMD`),
so the two are guaranteed consistent.

## How the samples were built

1. **Background — real, downloaded, never repeated.** Each sample is a distinct **1000 Genomes DRAGEN
   v4.4.7** exome, streamed from the public S3 bucket `1000genomes-dragen-v4-4-7`
   (`data/individuals/hg38-alt_masked.cnv.graph.hla.methyl_cg.rna-11-r5.0-2/<sample>/<sample>.hard-filtered.vcf.gz`),
   normalized (`bcftools norm -m -any`) and **subset to the MANE / GENCODE exome BED**
   (`data/gnomad/exome_hg38.bed`, ~100k variants each). Every case uses a **different** sample, chosen
   across diverse populations; no background is reused.
2. **Plant — one exact variant from a real patient case.** The causative variant is the **exact
   `chrom:pos:ref:alt`** recorded in a **GA4GH Phenopacket Store 0.1.27** case, together with that
   case's **HPO terms**, gene, consequence, and disease. When the coordinate is present in ClinVar, its
   real `CLNSIG/CLNREVSTAT/CLNDN` are carried; when it is **not** in the ClinVar release, a **synthetic
   `CLNSIG`** is used (disclosed in the manifest — those plants have no independent ClinVar assertion).
3. **Faithful genotype.** Each plant carries the patient's **real zygosity** from the phenopacket:
   **40 compound-heterozygous** (both true alleles), **75 homozygous**, **85 single-allele** — so
   recessive cases are biallelic, not lone carriers.
4. **Consequence spread.** Deliberately stratified across missense / stop-gained / frameshift /
   in-frame / start-loss (the expansion over-samples the missense/in-frame classes).

## The realistic (tell-free) transform — and the adjustments it needed

The default plant is trivially identifiable (`SPIKED=1` + `GENE/CSQ/CLN*` INFO + a minimal
`GT:DP:GQ:AD` FORMAT). The **realistic** build removes every marker: it **borrows a real background
call's full DRAGEN INFO/FORMAT** of the same zygosity, relocates it to the plant's coordinate/alleles,
and strips all markers. Truth is tracked **externally by coordinate** (the manifest / sidecars), never
in the VCF.

Two adjustments were required and are disclosed honestly:

- **QC-passing templates.** A borrowed real call can have a low genotype quality; the plant would then
  inherit it and be dropped at a caller/engine's QC. The template is required to pass QC comfortably
  (**GQ ≥ 30, DP ≥ 25, balanced het**). (Caught when a first build silently lost 17 cases.)
- **De-circularization is a feature.** Because the tell-free VCF carries **no** synthetic `CSQ`/`CLNSIG`,
  a tool must work from the **real sequence context**, not from the spike's hints. A handful of plants
  whose *marked* classification leaned on a synthetic tag (e.g. a `CSQ=missense` that the real annotator
  scores as a low-impact splice-region variant, or a `CLNSIG` for a variant absent from real ClinVar)
  land on their honest tier in the realistic build — a **cleaner, less circular** benchmark.

## Using it

- **Annotator validation:** annotate `realistic/*.vcf.gz`, then check each planted coordinate (from the
  manifest) got the expected gene / consequence / HGVS.
- **Caller validation:** the plant is a genuine-looking call; check it is recovered.
- **Report-generator validation:** run your pipeline on a sample + its `sidecars/SYN-NNN.hpo.txt`; the
  expected finding is the planted gene (see the manifest). Score blind on the `realistic/` set.

Reference: analysed with the [vcf2report](https://github.com/gbbarra/vcf2report) ACMG engine, the
planted variant reaches the diagnostic (primary) finding in **178/200** cases; the rest are honest
limitations (non-coding-RNA plants, HPO-unlinked genes, sub-threshold phenotype, missense held at VUS
without corroboration) — documented, not hidden.

## Honest limitations

- Synthetic: real backgrounds, **inserted** causative variant — tests classification/annotation, not a
  real diagnostic yield.
- Some plants carry a **synthetic ClinVar label** (coordinate absent from the ClinVar release) — flagged
  in the manifest.
- GRCh38 only.

## Sources & citations

- **1000 Genomes / IGSR**, DRAGEN v4.4.7 re-analysis (`1000genomes-dragen-v4-4-7`, AWS Open Data).
- **GA4GH Phenopacket Store** (Danis et al.) — the causative variants + HPO.
- **ClinVar** (NCBI, public domain) · **gnomAD v4.1** · **MANE / GENCODE** · **Human Phenotype Ontology (HPO)**.
- Annotation: **SnpEff 5.4c** (`GRCh38.mane.1.5.refseq`). Missense calibration in the reference engine
  uses **AlphaMissense** (CC BY-NC-SA — **not redistributed here**).

Tooling: MIT. Data honors each upstream source's license.
