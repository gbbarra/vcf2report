"""VCF parsing.

A dependency-free reader is the default so the pipeline runs anywhere (tests,
headless demo). If ``cyvcf2`` is installed it is used for speed/robustness on
real exomes. Either way we split multi-allelic records into single-allele
:class:`Variant`s and share one zygosity helper so the two paths cannot diverge.
NOTE: this does not left-align/trim indels — full normalization (``bcftools
norm -f REF``) is a production pre-step so keys match ClinVar/gnomAD.

Annotation carried in INFO is read from these keys when present:
``GENE``, ``CSQ`` (consequence), ``HGVSC``, ``HGVSP``. Real exomes annotated
with VEP/SnpEff can be mapped here; the bundled sample uses these plain keys.
"""
from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator, Optional

from ..models import Variant
from . import annparse

_MISSING = {".", "-1"}


def _open(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def detect_build(header_lines: list[str]) -> str | None:
    """Detect genome build from the VCF header using explicit, anchored tokens.

    Returns None when the build is not clearly declared so the pipeline emits its
    "build not declared" warning — a loose guess (e.g. any header containing the
    substring "38") could silently mislabel a GRCh37 VCF and defeat the guard.
    """
    blob = "\n".join(header_lines).lower()
    if "grch38" in blob or "hg38" in blob:
        return "GRCh38"
    if "grch37" in blob or "hg19" in blob or "b37" in blob or "g1k_v37" in blob:
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


def _alleles(gt: str) -> list[str]:
    return gt.replace("|", "/").split("/")


def zygosity(alleles: list[str], alt_num: int) -> Optional[str]:
    """Zygosity of the sample for a specific ALT allele number (1-based).

    Returns None for a no-call, a non-carrier (hom-ref, or carrying a *different*
    ALT), so callers can drop non-carriers. Shared by the pure and cyvcf2 paths.
    """
    alleles = [str(a) for a in alleles]
    if len(alleles) == 1:                      # hemizygous (e.g. chrX/Y male)
        a = alleles[0]
        if a in _MISSING:
            return None
        return "hemi" if a == str(alt_num) else None
    if any(a in _MISSING for a in alleles):
        return None                            # no-call -> unknown, never "hom"
    count = sum(1 for a in alleles if a == str(alt_num))
    if count == 0:
        return None                            # not a carrier of THIS allele
    return "hom" if count == len(alleles) else "het"  # 1/2 compound -> het


def _sample_metrics(fmt: str, sample: str, alt_index: int) -> dict:
    """Per-sample + per-ALT metrics. alt_index is 0-based into ALT."""
    d = dict(zip(fmt.split(":"), sample.split(":")))
    out: dict = {"zygosity": zygosity(_alleles(d.get("GT", "./.")), alt_index + 1)}
    if d.get("DP", ".").isdigit():
        out["depth"] = int(d["DP"])
    gq = d.get("GQ", ".")
    if gq not in _MISSING:
        try:
            out["gq"] = int(float(gq))
        except ValueError:
            pass
    ad = d.get("AD")
    if ad and "," in ad:
        try:
            parts = [int(x) for x in ad.split(",")]
            total = sum(parts)
            # AD is [ref, alt1, alt2, ...]; use THIS allele's depth.
            if total > 0 and len(parts) > alt_index + 1:
                out["allele_balance"] = round(parts[alt_index + 1] / total, 3)
        except ValueError:
            pass
    return out


def parse_vcf(path: str | Path) -> tuple[list[Variant], str | None, list[str]]:
    """Parse a VCF into (variants, detected_build, header_lines).

    Multi-allelic records are split into one :class:`Variant` per ALT allele,
    each carrying that allele's own zygosity and allele balance. Uses cyvcf2 if
    available, else the pure reader.
    """
    path = Path(path)
    try:  # pragma: no cover - exercised only when cyvcf2 present
        import cyvcf2  # type: ignore  # noqa: F401

        return _parse_cyvcf2(path)
    except Exception:
        return _parse_pure(path)


def _parse_pure(path: Path) -> tuple[list[Variant], str | None, list[str]]:
    variants: list[Variant] = []
    header: list[str] = []
    csq_format = None
    csq_done = False
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line)
                continue
            if not csq_done:  # header is complete once data starts
                csq_format = annparse.parse_csq_format(header)
                csq_done = True
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, filt, info = cols[:8]
            if not pos.isdigit() or not ref or not alt:
                continue  # malformed record — skip rather than crash
            info_d = _parse_info(info)
            fmt = cols[8] if len(cols) > 8 else ""
            sample = cols[9] if len(cols) > 9 else ""
            for i, alt_allele in enumerate(alt.split(",")):  # split multiallelics
                metrics = _sample_metrics(fmt, sample, i) if fmt and sample else {}
                ann = annparse.extract(info_d, alt_allele, csq_format, ref) or {}
                variants.append(Variant(
                    chrom=chrom, pos=int(pos), ref=ref, alt=alt_allele,
                    gene=ann.get("gene"),
                    hgvs_c=ann.get("hgvs_c"),
                    hgvs_p=ann.get("hgvs_p"),
                    consequence=ann.get("consequence"),
                    filter_status=filt,
                    zygosity=metrics.get("zygosity"),
                    depth=metrics.get("depth"),
                    gq=metrics.get("gq"),
                    allele_balance=metrics.get("allele_balance"),
                    info=info_d,
                    alt_index=i,
                ))
    return variants, detect_build(header), header


def _parse_cyvcf2(path: Path) -> tuple[list[Variant], str | None, list[str]]:  # pragma: no cover
    from cyvcf2 import VCF  # type: ignore

    vcf = VCF(str(path))
    header = [str(h) for h in vcf.raw_header.splitlines()]
    csq_format = annparse.parse_csq_format(header)
    variants: list[Variant] = []
    for rec in vcf:
        gts = rec.genotypes[0][:2] if rec.genotypes else None
        dp_arr = rec.format("DP")
        gq_arr = rec.format("GQ")
        ad_arr = rec.format("AD")
        depth = int(dp_arr[0][0]) if dp_arr is not None else None
        gq = None
        if gq_arr is not None:
            try:
                gq = int(gq_arr[0][0])
            except (ValueError, TypeError):
                gq = None
        for i, alt_allele in enumerate(rec.ALT):
            zyg = zygosity([str(a) for a in gts], i + 1) if gts else None
            allele_balance = None
            if ad_arr is not None:
                try:
                    ad = list(ad_arr[0])
                    total = sum(x for x in ad if x is not None and x >= 0)
                    if total > 0 and len(ad) > i + 1:
                        allele_balance = round(ad[i + 1] / total, 3)
                except (TypeError, ValueError, IndexError):
                    allele_balance = None
            info = {k: str(v) for k, v in dict(rec.INFO).items()}
            ann = annparse.extract(info, str(alt_allele), csq_format, rec.REF) or {}
            variants.append(Variant(
                chrom=rec.CHROM, pos=rec.POS, ref=rec.REF, alt=alt_allele,
                gene=ann.get("gene"), hgvs_c=ann.get("hgvs_c"),
                hgvs_p=ann.get("hgvs_p"), consequence=ann.get("consequence"),
                filter_status=rec.FILTER or "PASS", zygosity=zyg,
                depth=depth, gq=gq, allele_balance=allele_balance, info=info,
                alt_index=i,
            ))
    return variants, detect_build(header), header


def iter_variants(path: str | Path) -> Iterator[Variant]:
    variants, _, _ = parse_vcf(path)
    yield from variants
