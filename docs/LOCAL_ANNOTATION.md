# Native offline gnomAD + GRCh37→38 liftover (built in)

The pipeline's population-frequency step needs GRCh38 gnomAD. The **full** gnomAD v4.1
exomes+genomes sites set is **~150–200 GB** — usually too big for a laptop. So instead
of downloading all of it, vcf2report builds a **reduced tabix table** carrying only the
fields the ACMG engine cites (grpmax AF, AC/AN, homozygotes, the ClinGen filtering AF
`faf95`), one row per gnomAD variant. This is the piece that turns a real exome from
"11k VUS" into a short candidate list — the rarity filter can finally drop the common
variants and BA1/BS1/BS2 can fire. Only `pysam` is needed (already a dep).

## End-to-end on a real exome

```bash
# 0. (only if your VCF is GRCh37/hg19) lift it to GRCh38 first
VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/liftover_to_grch38.py exome.hg19.vcf exome.hg38.vcf

# 1. build a per-VCF local gnomAD table (tiny — one lookup per site, reuses the
#    same grpmax/faf95 reduction the live remote path uses, so local == remote)
VCF2REPORT_ALLOW_NETWORK=1 python3 scripts/build_gnomad_local.py --from-vcf exome.hg38.vcf
#    -> data/gnomad/gnomad_freq.local.tsv.gz (+ .tbi + .meta)

# 2. run fully offline — the report now narrows and calls Benign/VUS/Pathogenic
python3 scripts/run_headless.py exome.hg38.vcf --hpo hpo.txt --out out/
```

`gnomad.lookup` prefers the local table automatically when it exists; with no table
present, behaviour is unchanged (remote-if-network, else the bundled demo slice).
Point the table elsewhere (e.g. an external disk) with `VCF2REPORT_GNOMAD_TABIX=/path`.

## Build modes

| Mode | Command | Size / cost | Reusable for |
|---|---|---|---|
| **per-VCF** (default) | `--from-vcf V.vcf` | tiny; one lookup per site | that VCF's sites |
| **panel** | `--bed panel.bed` | small; the region's variants | any sample, those regions |
| **full** | `--full [--src DIR]` | **~150–200 GB streamed** (or local per-chrom sites files); genome-wide table | everything |

`--from-vcf` and `--bed` are the disk-friendly common cases. `--full` is for a machine
with room (use `--src` to read pre-downloaded per-chromosome sites files, no network).

## Safety model — never a false absence

A wrong "absent from gnomAD" (af 0.0) would make **PM2 fire** and **BA1/BS1 not fire**,
inflating pathogenicity. The table is built to make that impossible, and this was
verified by a 23-agent adversarial review:

- A **partial** table (per-VCF / panel) returns a value only on an **exact** variant
  match; any miss returns *unknown* and the caller falls back — it never asserts an
  absence. Genuinely-absent input sites still resolve via their explicit `af 0.0` row.
- A **full** table asserts absence only for a contig it **actually finished streaming**
  (recorded in the `.meta` sidecar's `contigs`). A chromosome whose stream failed, or a
  contig it structurally never covers (chrM/MT, alt/decoy), returns *unknown* → fallback,
  not a fabricated 0.0. An incomplete `--full` build warns loudly.
- The liftover drops ambiguous multi-target lifts and flags that it does not
  re-validate REF against the GRCh38 reference (normalize with `bcftools norm -f` if you
  have the FASTA).

---

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
