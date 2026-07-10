# Local annotation with SnpEff (run on your machine)

Consequence/HGVS + population/clinical annotation is a **local** step (no per-variant
network calls). It can't run in the sandbox (SnpEff's jar and GRCh38 database live on
hosts the build proxy blocks), but it's fast and offline on your machine. This turns a
raw / synthetic exome VCF into the fully-annotated VCF `vcf2report` reads directly.

## 1. Tools (install once — all free/MIT)

| Tool | What for | Install |
|---|---|---|
| **Java 11+** | runs SnpEff | `apt install default-jre` / `brew install openjdk` |
| **SnpEff** | consequence + HGVS (`ANN=`) | `curl -L https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip -o snpEff.zip && unzip snpEff.zip` (or `conda install -c bioconda snpeff`) |
| **bcftools + htslib** | normalize, bgzip, tabix | `conda install -c bioconda bcftools htslib` |
| **vcfanno** | gnomAD/ClinVar/REVEL from local files | `conda install -c bioconda vcfanno` |

Easiest: `conda create -n vcf2report -c bioconda -c conda-forge snpeff bcftools htslib vcfanno openjdk`

## 2. SnpEff GRCh38 database (~1–2 GB, one-time)

```bash
snpEff download GRCh38.105          # Ensembl-based; chromosomes are "1".."22","X"
```

> **Chromosome-naming gotcha (important for the DRAGEN synthetic cases).**
> DRAGEN VCFs are `chr`-prefixed (`chr1`), but the Ensembl SnpEff DB uses `1`. If they
> don't match, SnpEff annotates **nothing**. Pick ONE naming and make every file agree:
> - strip `chr` from the VCF before SnpEff:
>   `bcftools annotate --rename-chrs <(for c in $(seq 1 22) X Y M; do echo "chr$c $c"; done) in.vcf.gz -Oz -o in.nochr.vcf.gz`
> - use the Ensembl DB (`GRCh38.105`), and use a **no-`chr`** ClinVar (NCBI ships `1`-style)
>   and a **no-`chr`** gnomAD (or rename gnomAD, which is `chr`-style, the same way).

## 3. Annotation data files for vcfanno (GRCh38, bgzipped + tabixed)

Referenced by `scripts/vcfanno.conf.toml` — point each `file=` at your local path:

| File | Source | Notes |
|---|---|---|
| gnomAD v4.1 sites VCF | `gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/` (or genomes) | big; `chr`-style. Subset to your exome BED to shrink. |
| ClinVar GRCh38 VCF | `ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` | `1`-style, tabixed. (Sandbox-reachable mirror: the Broad Funcotator copy — uncompressed; `bgzip`+`tabix` it.) |
| ABraOM (optional, 🇧🇷) | `http://abraom.ib.usp.br/` | convert release to bgzipped VCF → the Brazilian-frequency differentiator |
| REVEL/CADD (optional) | dbNSFP / REVEL / CADD | in-silico predictors for PP3/BP4 |

## 4. Run it

```bash
# raw/synthetic VCF (bgzipped) + GRCh38 FASTA (+ .fai) -> annotated VCF
scripts/annotate_vcf.sh SYN-001.synthetic.vcf.gz GRCh38.fa SYN-001.annotated.vcf.gz
#   [1/3] bcftools norm   (split multiallelics + left-align)
#   [2/3] snpEff          (adds INFO/ANN: consequence + HGVS)   <-- SnpEff here
#   [3/3] vcfanno         (adds gnomAD_AF, CLNSIG, REVEL from local files)

# then the report — now the WHOLE background is annotated, so the candidate list is rich
python scripts/run_headless.py SYN-001.annotated.vcf.gz --hpo synthetic_exomes/SYN-001.hpo.txt
```

`SNPEFF_DB` (default `GRCh38.105`), `VCFANNO_CONF`, and `THREADS` are overridable env vars.

## Why not in the sandbox?

`snpeff.blob.core.windows.net` (jar + prebuilt DB) and every GTF mirror needed to
*build* a DB are blocked by the build proxy, and per-variant remote annotation is far
too slow. So annotation is deliberately a local step. The in-sandbox demo instead
spikes real ClinVar variants (which carry their own consequence) so the reports still
show correct Pathogenic/Likely-Pathogenic calls without the full background annotation.
