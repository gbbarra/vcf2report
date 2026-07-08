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
| 003 | NA18939 (EAS) | SCN2A | STK11 |
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

### In-sandbox alternative (no AWS CLI / bcftools / FASTA)

`scripts/build_dragen_exome.py` builds the real background exome purely over remote
tabix — it reads a **DRAGEN 4.4.7** per-sample VCF (bucket `1000genomes-dragen-v4-4-7`)
by HTTPS range queries and keeps PASS variants inside the IDT BED, no 440 MB
download and no bcftools. A full NA12878 exome (~24.8k real variants) builds in
~30 s. It was used to generate `docs/EXAMPLE_REPORT_SYN-001.md` end-to-end:

```bash
curl -fsSL "$IDT_BED_URL" -o scratch/idt_exome_v2.bed      # verified mirror (see make_synthetic_exomes.sh)
VCF2REPORT_ALLOW_NETWORK=1 python scripts/build_dragen_exome.py --sample NA12878 \
    --bed scratch/idt_exome_v2.bed --out data/real/NA12878_exome.vcf
# ClinVar subset for the 10 target genes (stream the reachable Broad mirror):
curl -s "$CLINVAR_MIRROR" | awk '/^#/ || /GENEINFO=(SCN1A|KCNQ2|SCN2A|STXBP1|SLC2A1|RB1|APC|STK11|WT1|FBN1):/' > scratch/clinvar_targets.vcf
python scripts/spike_pathogenic.py --exome data/real/NA12878_exome.vcf \
    --clinvar scratch/clinvar_targets.vcf --targets data/synthetic/SYN-001.targets.tsv \
    --sample-id SYN-001 --out data/real/SYN-001.synthetic.vcf
python -m vcf2report.cli data/real/SYN-001.synthetic.vcf --hpo data/synthetic/SYN-001.hpo.txt --out data/out
```

Note: without SnpEff in-sandbox the background variants carry no consequence and are
filtered at the impact stage, so the report surfaces the spiked P/LP findings amid
the real ~24.8k background; on a machine with the full annotation step the whole
background is annotated and the candidate list is richer.

### Verified results — all five cases built + reported end-to-end in-sandbox

Real DRAGEN 4.4.7 exome per sample, 2 real ClinVar pathogenic variants spiked,
gnomAD absence confirmed live, HPO-driven primary/secondary split. Reports in
`docs/example_reports/`.

| Case | Sample (super-pop) | Real variants | Primary (DEE) | Secondary (ACMG SF) |
|---|---|---|---|---|
| SYN-001 | NA12878 (EUR) | 24,801 | SCN1A → Pathogenic | RB1 → Pathogenic |
| SYN-002 | NA19240 (AFR) | 30,195 | KCNQ2 → Pathogenic | APC → Likely Pathogenic |
| SYN-003 | NA18939 (EAS) | 24,707 | SCN2A → Pathogenic | STK11 → Pathogenic |
| SYN-004 | NA20845 (SAS) | 25,319 | STXBP1 → Pathogenic | WT1 → Pathogenic |
| SYN-005 | HG01112 (AMR) | 24,387 | SLC2A1 → Pathogenic | FBN1 → Pathogenic |

## Notes
- **Chromosome naming** is harmonized to the exome's style (DRAGEN hg38 = `chr2`;
  ClinVar = `2`) by the spike script; keep chr-style references through annotation.
- Subsetting DRAGEN WGS to exome regions is a good demo proxy; the coverage profile
  differs from a real capture, which doesn't affect classification.
- **Not clinical data.** Synthetic, de-identified, marked demonstration-only in
  every output header.
