# SYN validation cohort (SYN-001..100)

A reproducible validation corpus of **100 real, distinct exomes** — each a different **1000G DRAGEN
v4.4.7** sample subset to the **vendor-neutral MANE/GENCODE exome BED** (`data/gnomad/exome_hg38.bed`,
the same region the engine covers) — with **one distinct pathogenic ClinVar variant spiked in**,
whose gene comes from a real **GA4GH phenopacket case** so it carries that case's **HPO** terms.

- **No repetition:** 100 distinct backgrounds, 100 distinct planted variants/genes. Different variants
  in different genomes → captures the planted signal across diverse backgrounds, simulating real cases.
- **SYN-001 = NA12878** (GIAB). The rest span the 1000G samples in the bucket.
- Variants: 25 stop-gain · 25 frameshift · 35 missense · 10 in-frame indel · 5 start-loss (all
  SNV/indel ≤50 bp — no SV/CNV). Source pool: 5,335 phenopacket variants-with-HPO (v0.1.27).

## Files

- **`cohort.tsv`** — the config (committed): `syn_id sample gene chrom pos ref alt consequence disease hpo`.
- `SYN-00N.synthetic.vcf.gz` (+ `.tbi`), `SYN-00N.hpo.txt`, `truth.tsv` — **generated** (git-ignored;
  rebuild anytime from `cohort.tsv`).

## Build (on a machine with AWS S3 access — not the sandbox)

```bash
REF_FASTA=/path/GRCh38.fa CLINVAR_VCF=/path/clinvar_GRCh38.vcf.gz \
  bash scripts/make_syn_cohort.sh
```
Needs: `awscli`, `bcftools`, `bgzip`, `tabix`, `python3`, a GRCh38 FASTA (+`.fai`), a ClinVar GRCh38 VCF.
Each case ≈ download + MANE subset + spike (~30 s–1 min); resumable (skips cases already built).

## Validate

```bash
bash scripts/validate_cohort.sh      # runs each case, checks the planted gene is surfaced vs truth.tsv
```
