#!/usr/bin/env python3
"""Provenance + regenerator for the bundled demo VCF.

Production recipe (documented; needs network + tools):
  1. Download a public, de-identified exome — e.g. Genome in a Bottle NA12878
     (HG001), GRCh38, from the GIAB FTP.
  2. Subset to a handful of genes with clean HPO annotations to keep it tiny.
  3. Spike 2-4 known ClinVar pathogenic variants that match one chosen
     rare-disease phenotype so the demo reliably shows P/LP calls.
  4. Strip all sample identifiers; keep only an opaque ID (DEMO-001).

For this repo we ship a small, fully SYNTHETIC, de-identified VCF (no real
patient data). Running this script rewrites data/sample/sample_exome.vcf
deterministically so the demo is reproducible without any download.

    python scripts/build_sample_vcf.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from vcf2report import config  # noqa: E402

HEADER = """##fileformat=VCFv4.2
##fileDate=2024-11-01
##source=vcf2report_demo
##reference=GRCh38
##contig=<ID=1>
##contig=<ID=2>
##contig=<ID=3>
##contig=<ID=7>
##contig=<ID=13>
##contig=<ID=11>
##contig=<ID=19>
##contig=<ID=20>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CSQ,Number=1,Type=String,Description="Molecular consequence">
##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS coding">
##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVS protein">
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##comment=SYNTHETIC, de-identified demonstration variants. Not real patient data. Not for clinical use.
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO-001
"""

# (chrom, pos, ref, alt, gene, csq, hgvsc, hgvsp, filt, sample)
ROWS = [
    ("2", 166003360, "C", "T", "SCN1A", "stop_gained", "c.1834C>T", "p.Arg612Ter", "PASS", "0/1:45:99:22,23"),
    ("20", 63446204, "G", "A", "KCNQ2", "missense_variant", "c.637C>T", "p.Arg213Trp", "PASS", "0/1:38:99:19,19"),
    ("19", 13318000, "T", "C", "CACNA1A", "missense_variant", "c.100T>C", "p.Ser34Pro", "PASS", "0/1:40:99:20,20"),
    # Incidental Likely Pathogenic: LoF in a LoF-intolerant gene unrelated to the
    # seizure indication (no phenotype match -> no PP4, no ClinVar -> no PP5) so it
    # lands at LP (PVS1 + PM2), not over-called to Pathogenic.
    ("11", 31806340, "G", "A", "PAX6", "stop_gained", "c.202C>T", "p.Arg68Ter", "PASS", "0/1:46:99:23,23"),
    # Incidental finding on an ACMG SF v3.2 gene (RB1 / hereditary retinoblastoma):
    # a reportable, actionable secondary finding, cleanly unrelated to the seizure
    # indication (RB1 shares no HPO term with the patient's seizure phenotype) and
    # LoF-intolerant, so PVS1 + PM2 -> Likely Pathogenic.
    ("13", 48367226, "C", "T", "RB1", "stop_gained", "c.958C>T", "p.Arg320Ter", "PASS", "0/1:44:99:22,22"),
    ("1", 228208000, "G", "A", "OBSCN", "missense_variant", "c.298G>A", "p.Val100Ile", "PASS", "0/1:50:99:26,24"),
    ("2", 178562809, "G", "A", "TTN", "missense_variant", "c.1000G>A", "p.Val334Ile", "PASS", "0/1:60:99:30,30"),
    ("1", 11790000, "A", "G", "MTOR", "intron_variant", "", "", "PASS", "0/1:42:99:21,21"),
    ("3", 37000000, "C", "T", "MLH1", "synonymous_variant", "c.655C>T", "p.=", "PASS", "0/1:48:99:24,24"),
    ("7", 117559590, "G", "T", "CFTR", "missense_variant", "c.254G>T", "p.Gly85Val", "PASS", "0/1:6:30:3,3"),
    ("11", 5227000, "A", "T", "HBB", "missense_variant", "c.20A>T", "p.Glu7Val", "LowQual", "0/1:35:45:18,17"),
]


def main() -> int:
    lines = [HEADER.rstrip("\n")]
    for chrom, pos, ref, alt, gene, csq, hgvsc, hgvsp, filt, sample in ROWS:
        info = f"GENE={gene};CSQ={csq}"
        if hgvsc:
            info += f";HGVSC={hgvsc}"
        if hgvsp:
            info += f";HGVSP={hgvsp}"
        lines.append("\t".join([chrom, str(pos), ".", ref, alt, "800", filt,
                                info, "GT:DP:GQ:AD", sample]))
    config.SAMPLE_VCF.parent.mkdir(parents=True, exist_ok=True)
    config.SAMPLE_VCF.write_text("\n".join(lines) + "\n")
    print(f"Wrote {len(ROWS)} variants to {config.SAMPLE_VCF}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
