"""VCF parsing.

A dependency-free reader is the default so the pipeline runs anywhere (tests,
headless demo). If ``cyvcf2`` is installed it is used for speed/robustness on
real exomes. Either way we emit normalized single-allele :class:`Variant`s.

Annotation carried in INFO is read from these keys when present:
``GENE``, ``CSQ`` (consequence), ``HGVSC``, ``HGVSP``. Real exomes annotated
with VEP/SnpEff can be mapped here; the bundled sample uses these plain keys.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator

from ..models import Variant


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def detect_build(header_lines: list[str]) -> str | None:
    """Best-effort genome-build detection from the VCF header."""
    blob = "\n".join(header_lines).lower()
    if "grch38" in blob or "hg38" in blob or "38" in blob and "reference" in blob:
        return "GRCh38"
    if "grch37" in blob or "hg19" in blob:
        return "GRCh37"
    return None


def _parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if info in (".", ""):
        return out
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = "true"
    return out


def _zygosity(gt: str) -> str | None:
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) == 1:
        return "hemi"
    if "." in alleles:
        return None
    nonref = [a for a in alleles if a not in ("0",)]
    if not nonref:
        return None
    return "hom" if len(set(alleles)) == 1 else "het"


def _sample_metrics(fmt: str, sample: str) -> dict:
    keys = fmt.split(":")
    vals = sample.split(":")
    d = dict(zip(keys, vals))
    out: dict = {"zygosity": _zygosity(d.get("GT", "./."))}
    if d.get("DP", ".").isdigit():
        out["depth"] = int(d["DP"])
    if d.get("GQ", ".").replace(".", "", 1).isdigit() and d.get("GQ") != ".":
        try:
            out["gq"] = int(float(d["GQ"]))
        except ValueError:
            pass
    ad = d.get("AD")
    if ad and "," in ad:
        try:
            parts = [int(x) for x in ad.split(",")]
            total = sum(parts)
            if total > 0 and len(parts) >= 2:
                out["allele_balance"] = round(parts[1] / total, 3)
        except ValueError:
            pass
    return out


def parse_vcf(path: str | Path) -> tuple[list[Variant], str | None, list[str]]:
    """Parse a VCF into (variants, detected_build, header_lines).

    Multi-allelic records are split into one :class:`Variant` per ALT allele
    (basic normalization). Uses cyvcf2 if available, else the pure reader.
    """
    path = Path(path)
    try:  # pragma: no cover - exercised only when cyvcf2 present
        from cyvcf2 import VCF  # type: ignore

        return _parse_cyvcf2(path)
    except Exception:
        return _parse_pure(path)


def _parse_pure(path: Path) -> tuple[list[Variant], str | None, list[str]]:
    variants: list[Variant] = []
    header: list[str] = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, filt, info = cols[:8]
            info_d = _parse_info(info)
            fmt = cols[8] if len(cols) > 8 else ""
            sample = cols[9] if len(cols) > 9 else ""
            metrics = _sample_metrics(fmt, sample) if fmt and sample else {}
            for alt_allele in alt.split(","):  # split multiallelics
                variants.append(Variant(
                    chrom=chrom, pos=int(pos), ref=ref, alt=alt_allele,
                    gene=info_d.get("GENE"),
                    hgvs_c=info_d.get("HGVSC"),
                    hgvs_p=info_d.get("HGVSP"),
                    consequence=info_d.get("CSQ"),
                    filter_status=filt,
                    zygosity=metrics.get("zygosity"),
                    depth=metrics.get("depth"),
                    gq=metrics.get("gq"),
                    allele_balance=metrics.get("allele_balance"),
                ))
    return variants, detect_build(header), header


def _parse_cyvcf2(path: Path) -> tuple[list[Variant], str | None, list[str]]:  # pragma: no cover
    from cyvcf2 import VCF  # type: ignore

    vcf = VCF(str(path))
    header = [str(h) for h in vcf.raw_header.splitlines()]
    variants: list[Variant] = []
    for rec in vcf:
        for i, alt_allele in enumerate(rec.ALT):
            info = dict(rec.INFO)
            gt = rec.genotypes[0][:2] if rec.genotypes else None
            zyg = None
            if gt:
                zyg = "hom" if gt[0] == gt[1] and gt[0] != 0 else ("het" if 0 in gt else "hom")
            variants.append(Variant(
                chrom=rec.CHROM, pos=rec.POS, ref=rec.REF, alt=alt_allele,
                gene=info.get("GENE"), hgvs_c=info.get("HGVSC"),
                hgvs_p=info.get("HGVSP"), consequence=info.get("CSQ"),
                filter_status=rec.FILTER or "PASS", zygosity=zyg,
                depth=rec.format("DP")[0][0] if rec.format("DP") is not None else None,
            ))
    return variants, detect_build(header), header


def iter_variants(path: str | Path) -> Iterator[Variant]:
    variants, _, _ = parse_vcf(path)
    yield from variants
