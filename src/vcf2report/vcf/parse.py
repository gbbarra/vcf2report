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
import os
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
            # Many callers omit FORMAT/DP but carry AD; sum(AD) is the site depth.
            if "depth" not in out and total > 0:
                out["depth"] = total
        except ValueError:
            pass
    return out


def _reportable_alt(alt: str) -> bool:
    """False for non-sequence ALTs that must never inherit another allele's
    annotation: spanning deletion '*', symbolic '<DEL>'/'<NON_REF>'/'<*>', '.'."""
    return bool(alt) and alt not in ("*", ".") and not alt.startswith("<")


def _sample_names(header: list[str]) -> list[str]:
    for line in header:
        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            return cols[9:] if len(cols) > 9 else []
    return []


def _resolve_sample_index(header: list[str], sample: str | None) -> int:
    names = _sample_names(header)
    if sample is None:
        return 0
    if sample in names:
        return names.index(sample)
    raise ValueError(f"sample {sample!r} not found in VCF (samples: {names})")


def parse_vcf(path: str | Path, sample: str | None = None
              ) -> tuple[list[Variant], str | None, list[str]]:
    """Parse a VCF into (variants, detected_build, header_lines).

    Multi-allelic records are split into one :class:`Variant` per ALT allele
    (non-sequence alleles like '*' / '<...>' are skipped), each carrying that
    allele's own zygosity and allele balance. For a multi-sample VCF, ``sample``
    selects the proband by name (defaults to the first column; the pipeline warns
    when a multi-sample VCF is parsed without an explicit selection).
    Uses cyvcf2 if available, else the pure reader.
    """
    path = Path(path)
    # cyvcf2 (htslib) is the fast path for real exomes; the pure reader is the
    # dependency-free default. Set VCF2REPORT_NO_CYVCF2=1 to force the pure reader.
    if os.environ.get("VCF2REPORT_NO_CYVCF2", "").strip().lower() not in {"1", "true", "yes"}:
        try:  # pragma: no cover - exercised only when cyvcf2 present
            import cyvcf2  # type: ignore  # noqa: F401

            return _parse_cyvcf2(path, sample)
        except Exception:
            pass
    return _parse_pure(path, sample)


def _parse_pure(path: Path, sample: str | None = None
                ) -> tuple[list[Variant], str | None, list[str]]:
    variants: list[Variant] = []
    header: list[str] = []
    csq_format = None
    sample_idx = 0
    setup_done = False
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")  # tolerate CRLF
            if line.startswith("#"):
                header.append(line)
                continue
            if not setup_done:  # header is complete once data starts
                csq_format = annparse.parse_csq_format(header)
                sample_idx = _resolve_sample_index(header, sample)
                setup_done = True
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, ref, alt, _qual, filt, info = cols[:8]
            if not pos.isdigit() or not ref or not alt:
                continue  # malformed record — skip rather than crash
            info_d = _parse_info(info)
            fmt = cols[8] if len(cols) > 8 else ""
            scol = 9 + sample_idx
            samp = cols[scol] if len(cols) > scol else ""
            alts = alt.split(",")
            for i, alt_allele in enumerate(alts):  # split multiallelics
                if not _reportable_alt(alt_allele):
                    continue  # '*' / symbolic ALT: never annotate/report
                metrics = _sample_metrics(fmt, samp, i) if fmt and samp else {}
                ann = annparse.extract(info_d, alt_allele, csq_format, ref, i, len(alts)) or {}
                variants.append(Variant(
                    chrom=chrom, pos=int(pos), ref=ref, alt=alt_allele,
                    gene=ann.get("gene"),
                    hgvs_c=ann.get("hgvs_c"),
                    hgvs_p=ann.get("hgvs_p"),
                    consequence=ann.get("consequence"),
                    exon=ann.get("exon"),
                    transcript=ann.get("transcript"),
                    filter_status=filt,
                    variant_id=_id if _id not in (".", "") else None,
                    n_alts=len(alts),
                    zygosity=metrics.get("zygosity"),
                    depth=metrics.get("depth"),
                    gq=metrics.get("gq"),
                    allele_balance=metrics.get("allele_balance"),
                    info=info_d,
                    alt_index=i,
                ))
    return variants, detect_build(header), header


def _cyvcf2_int(v):  # pragma: no cover
    """cyvcf2 encodes a missing per-sample INT as the INT32 min sentinel."""
    try:
        iv = int(v)
    except (ValueError, TypeError):
        return None
    return iv if iv >= 0 else None


def _parse_cyvcf2(path: Path, sample: str | None = None
                  ) -> tuple[list[Variant], str | None, list[str]]:  # pragma: no cover
    from cyvcf2 import VCF  # type: ignore

    vcf = VCF(str(path))
    header = [str(h) for h in vcf.raw_header.splitlines()]
    csq_format = annparse.parse_csq_format(header)
    s = 0
    if sample is not None:
        if sample not in list(vcf.samples):
            raise ValueError(f"sample {sample!r} not found (samples: {list(vcf.samples)})")
        s = list(vcf.samples).index(sample)
    variants: list[Variant] = []
    for rec in vcf:
        gts = rec.genotypes[s][:2] if rec.genotypes else None
        dp_arr = rec.format("DP")
        gq_arr = rec.format("GQ")
        ad_arr = rec.format("AD")
        depth = _cyvcf2_int(dp_arr[s][0]) if dp_arr is not None else None
        gq = _cyvcf2_int(gq_arr[s][0]) if gq_arr is not None else None
        alts = list(rec.ALT)
        for i, alt_allele in enumerate(alts):
            if not _reportable_alt(str(alt_allele)):
                continue
            zyg = zygosity([str(a) for a in gts], i + 1) if gts else None
            allele_balance = None
            depth_i = depth
            if ad_arr is not None:
                try:
                    ad = [x for x in list(ad_arr[s]) if x is not None and x >= 0]
                    total = sum(ad)
                    if total > 0 and len(ad) > i + 1:
                        allele_balance = round(ad[i + 1] / total, 3)
                    if depth_i is None and total > 0:  # DP absent -> sum(AD)
                        depth_i = total
                except (TypeError, ValueError, IndexError):
                    allele_balance = None
            info = {k: str(v) for k, v in dict(rec.INFO).items()}
            ann = annparse.extract(info, str(alt_allele), csq_format, rec.REF, i, len(alts)) or {}
            variants.append(Variant(
                chrom=rec.CHROM, pos=rec.POS, ref=rec.REF, alt=alt_allele,
                gene=ann.get("gene"), hgvs_c=ann.get("hgvs_c"),
                hgvs_p=ann.get("hgvs_p"), consequence=ann.get("consequence"),
                exon=ann.get("exon"), transcript=ann.get("transcript"),
                filter_status=rec.FILTER or "PASS", zygosity=zyg,
                variant_id=rec.ID if rec.ID not in (None, ".", "") else None,
                n_alts=len(alts),
                depth=depth_i, gq=gq, allele_balance=allele_balance, info=info,
                alt_index=i,
            ))
    return variants, detect_build(header), header


def iter_variants(path: str | Path) -> Iterator[Variant]:
    variants, _, _ = parse_vcf(path)
    yield from variants
