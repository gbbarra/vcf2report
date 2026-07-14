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

## Build

```bash
cd ~/vcf2report && bash scripts/make_syn_cohort.sh      # all 100 (resumable)
# smoke test first:  N=3 bash scripts/make_syn_cohort.sh
```
Self-contained — needs only `curl`, `bcftools`, `bgzip`, `tabix`, `python3`. The DRAGEN VCFs stream
from the public S3 bucket over HTTPS (no `aws` CLI); the ClinVar VCF is auto-downloaded once (~180 MB);
no reference FASTA. Each case ≈ 422 MB download + MANE subset + spike (~1–3 min); resumable (skips
cases already built). Point at an existing ClinVar VCF with `CLINVAR_VCF=/path/clinvar.vcf.gz` to skip
the download.

## Validate

```bash
bash scripts/validate_cohort.sh      # runs each case, checks the planted gene is surfaced vs truth.tsv
```
