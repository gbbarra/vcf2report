"""Build the local ClinVar store for offline, complete, random-access lookups.

Converts the official ClinVar GRCh38 VCF into a bgzipped, position-sorted TSV
(chrom pos ref alt significance review_status accession condition) indexed with tabix,
so vcf2report can look up any ClinVar variant offline without a network call and without
loading millions of rows into memory. Underscores in the text fields are normalised to
spaces (matching the E-utilities form that clinvar_stars/PP5 parse).

    python3 scripts/build_clinvar_local.py <clinvar_grch38.vcf.gz> [out.tsv.gz]

Default output: data/clinvar/clinvar_grch38.tsv.gz (+ .tbi). Get the source VCF from
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
"""
import gzip, os, re, subprocess, sys

SRC = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "clinvar", "clinvar_grch38.tsv.gz")
TMP = OUT[:-3] if OUT.endswith(".gz") else OUT + ".tsv"


def field(info, key):
    m = re.search(rf"(?:^|;){key}=([^;]+)", info)
    return m.group(1).replace("_", " ") if m else ""


rows, kept = [], 0
op = gzip.open if SRC.endswith(".gz") else open
with op(SRC, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 8:
            continue
        chrom, pos, vid, ref, alt, info = c[0], c[1], c[2], c[3], c[4], c[7]
        sig = field(info, "CLNSIG")
        if not sig or alt in (".", ""):
            continue
        rev = field(info, "CLNREVSTAT")
        cond = field(info, "CLNDN")[:120]
        rows.append((chrom, int(pos), ref.upper(), alt.upper(), sig, rev, vid, cond))
        kept += 1

rows.sort(key=lambda r: (r[0], r[1]))
os.makedirs(os.path.dirname(OUT), exist_ok=True)
with open(TMP, "w") as o:
    for r in rows:
        o.write("\t".join(map(str, r)) + "\n")
subprocess.run(["bgzip", "-f", TMP], check=True)
subprocess.run(["tabix", "-s", "1", "-b", "2", "-e", "2", "-f", OUT], check=True)
print(f"wrote {OUT} + .tbi | {kept:,} ClinVar variants (from {SRC})")
