# Example report — NA12878 (real whole exome, fully offline)

_Auto-generated from a real run; not hand-edited. This is a **concise excerpt** (funnel + performance + the non-VUS calls); the full run classifies every candidate. Reproduce with a local gnomAD Parquet store (`scripts/fetch_gnomad_parquet.sh` or `build_gnomad_parquet.py`) and:_

```bash
python3 -m vcf2report.cli NA12878_exome.hg38.vcf \
  --hpo HP:0001250,HP:0001263 --sample-id NA12878 --out out/
```

**Sample:** GIAB NA12878 / HG001 (NIST, public domain) · TruSeq exome · lifted to GRCh38.  
**Frequencies:** local gnomAD v4.1 DuckDB/Parquet store (29.6M variants), no network.

## Narrowing funnel

| stage | variants |
|---|---|
| total called | 28,565 |
| pass QC (DP/GQ/AB) | 23,773 |
| rare (gnomAD/ABraOM) | 4,656 |
| impactful / ClinVar P-LP | 1,564 |
| **candidates classified** | **1,564** |

Tiers: Benign 1, Likely Benign 60, Likely Pathogenic 3, Uncertain Significance (VUS) 1500

## Performance (this run)

| stage | seconds |
|---|---|
| parse | 0.812 |
| qc | 0.0351 |
| gnomad prime | 3.1119 |
| annotate | 3.6232 |
| filter | 0.0269 |
| alphamissense | 0.6848 |
| classify | 0.2255 |
| total | 8.5194 |
| **wall (incl. import)** | **8.6** |

~3352.9 variants/s, fully offline.

## Sequencing quality (proxy at variant sites)

Mean depth 64.0× (median 47.0×) · Ti/Tv 2.83 · het:hom 1.74 · assay exome / large-panel-scale

## Non-VUS calls (64)

| variant | gene | consequence | tier | criteria |
|---|---|---|---|---|
| 17-63583528-CG-C | DCAF7 | splice_acceptor_variant | Likely Pathogenic | PVS1, PM2 |
| 18-76379002-CG-C | ZNF516 | frameshift_variant | Likely Pathogenic | PVS1, PM2 |
| p.Trp31Arg | CDC27 | missense_variant | Likely Pathogenic | PM2, PP3 |
| p.Gly140Arg | PEX10 | missense_variant | Likely Benign | PP4, BS2, BP4 |
| p.Val289Ile | HDAC4 | missense_variant | Likely Benign | PP4, BS2, BP4 |
| p.Val625Met | VARS2 | missense_variant | Likely Benign | PP4, BS2, BP4 |
| p.Thr1205Ala | KCNT1 | missense_variant | Likely Benign | PP4, BS2, BP4 |
| p.Glu360Asp | GLIS3 | missense_variant | Likely Benign | PP4, BS2, BP4 |
| p.Ser571Leu | RBM28 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro219Ser | TPRN | missense_variant | Likely Benign | BS2, BP4 |
| p.Ser300Pro | CEP78 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg554Gln | DENND2A | missense_variant | Likely Benign | BS2, BP4 |
| p.Leu198Ile | OR4N4 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ile598Val | DGKG | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala454Thr | DOCK9 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro602Ser | VWA5A | missense_variant | Likely Benign | BS2, BP4 |
| p.Ile1447Thr | COL4A4 | missense_variant | Likely Benign | BS2, BP4 |
| p.Thr149Met | CSH2 | missense_variant | Likely Benign | BS2, BP4 |
| p.Phe608Tyr | SLC7A4 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro1330Ser | SSC5D | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg246Cys | TOP1MT | missense_variant | Likely Benign | BS2, BP4 |
| p.Ile170Thr | CD74 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ser358Cys | GADL1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Thr329Ala | HLA-C | missense_variant | Likely Benign | BS2, BP4 |
| p.Cys124Ser | PRAMEF6 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala330Val | HLA-C | missense_variant | Likely Benign | BS2, BP4 |
| p.Asp250Ala | POLQ | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala30Thr | PAPL | missense_variant | Likely Benign | BS2, BP4 |
| p.Met331Val | HLA-C | missense_variant | Likely Benign | BS2, BP4 |
| p.His62Tyr | DEFB125 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro685Leu | SALL3 | missense_variant | Likely Benign | BS2, BP4 |
| p.Val319Ala | HLA-C | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg244His | ABTB1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Leu324Val | PRSS53 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ile435Thr | GRIK1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Leu300Met | MARVELD2 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg1400Trp | MYH9 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro708Ser | ABCA3 | missense_variant | Likely Benign | BS2, BP4 |
| p.Asn100Ser | HOXC10 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro320Leu | SGIP1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg4Cys | KRT84 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala1003Thr | MYO10 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ser293Arg | PPAPDC2 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg113Gln | SAA4 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ile112Met | SH3D19 | missense_variant | Likely Benign | BS2, BP4 |
| p.Trp4Arg | ATP6AP1L | missense_variant | Likely Benign | BS2, BP4 |
| p.Thr361Ala | TXLNB | missense_variant | Likely Benign | BS2, BP4 |
| p.Met309Val | HLA-C | missense_variant | Likely Benign | BS2, BP4 |
| p.Leu424Phe | MAGEC1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Pro194Leu | FPGT | missense_variant | Likely Benign | BS2, BP4 |
| p.Asn1383Ser | ADAMTS12 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg319Cys | ARHGAP8 | missense_variant | Likely Benign | BS2, BP4 |
| p.Val503Ile | GPR113 | missense_variant | Likely Benign | BS2, BP4 |
| p.Thr208Met | FRMPD2 | missense_variant | Likely Benign | BS2, BP4 |
| p.Glu170Gln | RFPL3 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala42Val | STK19 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg34Gln | C20orf85 | missense_variant | Likely Benign | BS2, BP4 |
| p.Val210Met | FHAD1 | missense_variant | Likely Benign | BS2, BP4 |
| p.Arg23Lys | IL1F10 | missense_variant | Likely Benign | BS2, BP4 |
| p.Ala520Val | MYO3A | missense_variant | Likely Benign | BS2, BP4 |

_… 4 more not shown._

_The remaining 1,500 candidates are VUS — conservative by design when population frequency alone can't resolve them (no over-calling)._

