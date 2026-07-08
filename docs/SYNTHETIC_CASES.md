# Synthetic case exomes (real 1000G background + spiked ClinVar pathogenic)

The strongest demo input is a **real, diverse, exome-scale VCF that still reliably
shows Pathogenic/Likely-Pathogenic calls**. Raw 1000G samples are healthy, so on
their own they only yield benign/VUS — a weak demo. The fix: keep the **real
background** and plant **real ClinVar pathogenic variants** at their true GRCh38
coordinates. Everything is synthetic-but-real-sourced and de-identified.

```
1000G DRAGEN exome (real, healthy)  ─┐
                                     ├─ spike_pathogenic.py ─► SYN-00N.synthetic.vcf
ClinVar pathogenic @ true coords  ──┘        (de-identified, real coords/HGVS)
        │
        ├─ annotate_vcf.sh (SnpEff + vcfanno: consequence + gnomAD + ClinVar)
        ▼
   run_headless.py --hpo ...  ─►  auditable ACMG report (primary / secondary split)
```

## Two scripts

- **`scripts/spike_pathogenic.py`** — pulls the chosen genes' pathogenic records
  straight from the ClinVar VCF (LoF preferred, deterministic), plants them into a
  real exome VCF at their true coordinates, carries `CLNSIG`/`CLNREVSTAT`/`CLNDN`
  (so PP5 fires and the disease name reaches the report), and de-identifies the
  sample to `SYN-00N`.
- **`scripts/make_synthetic_exomes.sh`** — the full per-sample pipeline: download
  from the public 1000G DRAGEN S3 bucket → `bcftools norm` → subset to the IDT
  xGen Exome v2 targets → spike → sort/bgzip/tabix, plus the matching HPO file.

## The five validated cases

`data/synthetic/SYN-00N.{targets.tsv,hpo.txt}` — see `data/synthetic/README.md`.
Genes were chosen so each routes correctly against the **full real HPO** + **real
gnomAD constraint**: primary = LoF-intolerant DEE gene that matches the seizure
phenotype (→ Pathogenic, primary); secondary = LoF-intolerant ACMG SF gene with
zero seizure-HPO overlap (→ Pathogenic, secondary, not leaking into primary).

| SYN | 1000G (super-pop) | Primary (DEE) | Secondary (ACMG SF) |
|---|---|---|---|
| 001 | NA12878 (EUR) | SCN1A | RB1 |
| 002 | NA19240 (AFR) | KCNQ2 | APC |
| 003 | NA18939 (EAS) | SCN2A | SMAD4 |
| 004 | NA20845 (SAS) | STXBP1 | WT1 |
| 005 | HG01112 (AMR) | SLC2A1 | FBN1 |

## Run (on a machine that can reach S3 + ClinVar + a GRCh38 FASTA)

```bash
# 1. inputs: IDT xGen Exome v2 BED (hg38), ClinVar GRCh38 VCF, GRCh38.fa(+.fai), awscli
# 2. build all five synthetic exomes
scripts/make_synthetic_exomes.sh

# 3. annotate + report one case
scripts/annotate_vcf.sh synthetic_exomes/SYN-001.synthetic.vcf.gz GRCh38.fa SYN-001.annotated.vcf.gz
python scripts/run_headless.py SYN-001.annotated.vcf.gz --hpo synthetic_exomes/SYN-001.hpo.txt
```

The spike step alone is verifiable offline against a tiny mock ClinVar; the S3 /
FASTA / IDT-BED steps require the real inputs and are meant to run on your machine.

## Notes
- **Chromosome naming** is harmonized to the exome's style (DRAGEN hg38 = `chr2`;
  ClinVar = `2`) by the spike script; keep chr-style references through annotation.
- Subsetting DRAGEN WGS to exome regions is a good demo proxy; the coverage profile
  differs from a real capture, which doesn't affect classification.
- **Not clinical data.** Synthetic, de-identified, marked demonstration-only in
  every output header.
