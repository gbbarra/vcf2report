"""Live gnomAD frequencies via remote tabix over the public GCS bucket.

gnomAD publishes its per-chromosome sites VCFs (bgzipped + tabix-indexed) in the
``gcp-public-data--gnomad`` bucket. htslib can range-query them over HTTPS, so we
can pull the real grpmax (popmax) allele frequency for a single variant without
downloading the whole file and without the GraphQL API's rate limits.

This is the live path that actually works in locked-down networks where only
``storage.googleapis.com`` is reachable. It queries the current release
(v4.1 exomes + genomes) and reports the higher grpmax of the two callsets.

Open file handles are cached per chromosome (opening a remote VCF downloads the
header + tabix index, which is too costly to repeat per variant), so annotating a
whole exome reuses ~24 handles.
"""
from __future__ import annotations

import threading
from typing import Optional

from ..models import Variant

RELEASE = "4.1"
_BASE = ("https://storage.googleapis.com/gcp-public-data--gnomad/release/"
         f"{RELEASE}/vcf/{{kind}}/gnomad.{{kind}}.v{RELEASE}.sites.{{chrom}}.vcf.bgz")

_handles: dict[tuple[str, str], object] = {}
_failed: set[tuple[str, str]] = set()
_lock = threading.Lock()
_pysam = None
_pysam_tried = False


def _get_pysam():
    global _pysam, _pysam_tried
    if not _pysam_tried:
        _pysam_tried = True
        try:
            import pysam  # noqa
            _pysam = pysam
        except Exception:
            _pysam = None
    return _pysam


def _chrom(variant: Variant) -> str:
    c = str(variant.chrom)
    return c if c.startswith("chr") else f"chr{c}"


def _open(kind: str, chrom: str):
    """Return a cached VariantFile for (kind, chrom), or None if unavailable."""
    key = (kind, chrom)
    if key in _failed:
        return None
    h = _handles.get(key)
    if h is not None:
        return h
    with _lock:
        if key in _failed:
            return None
        h = _handles.get(key)
        if h is not None:
            return h
        pysam = _get_pysam()
        if pysam is None:
            _failed.add(key)
            return None
        url = _BASE.format(kind=kind, chrom=chrom)
        try:
            h = pysam.VariantFile(url)
        except Exception:
            _failed.add(key)
            return None
        _handles[key] = h
        return h


_POPMAX_EXCLUDE = {"asj", "fin", "oth", "remaining", "mid", "ami"}


def _best_from_record(rec) -> Optional[dict]:
    """grpmax (popmax) AF/AC/AN/hom for a matched single-ALT record.

    Falls back to the global AF/AC/AN when grpmax is undefined (e.g. the max lies
    in an excluded bottleneck group) so a real, non-absent variant is never
    reported as AF 0.
    """
    def a0(key):
        v = rec.info.get(key)
        if isinstance(v, (tuple, list)):
            return v[0] if v else None
        return v

    # Filtering AF (95% CI upper bound of the grpmax group) — the field ClinGen
    # recommends for BS1/BA1. gnomAD v4.1 publishes it site-wide as fafmax_faf95_max.
    def faf95():
        for k in ("fafmax_faf95_max", "faf95_grpmax", "faf95_max"):
            v = a0(k)
            if v is not None:
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
        return None

    pop = a0("grpmax")
    if pop is not None and pop.lower() not in _POPMAX_EXCLUDE:
        af = a0("AF_grpmax")
        if af is not None:
            return {"af": float(af), "ac": a0("AC_grpmax"), "an": a0("AN_grpmax"),
                    "hom": a0("nhomalt_grpmax"), "faf95": faf95(), "pop": pop}
    # fall back to global
    af = a0("AF")
    if af is None:
        return None
    return {"af": float(af), "ac": a0("AC"), "an": a0("AN"),
            "hom": a0("nhomalt"), "faf95": faf95(), "pop": None}


def query(variant: Variant) -> Optional[dict]:
    """Real gnomAD grpmax frequency for ``variant`` via remote tabix.

    Returns a dict ``{af, ac, an, hom, pop}`` (the higher grpmax of exomes vs
    genomes), a confirmed-absent dict (af 0) if the position is covered by gnomAD
    but the exact allele is not present, or ``None`` if the lookup could not run
    (pysam missing / all chromosome handles failed) so the caller can fall back.
    """
    pysam = _get_pysam()
    if pysam is None:
        return None
    chrom = _chrom(variant)
    pos = variant.pos
    opened = 0
    best: Optional[dict] = None
    for kind in ("exomes", "genomes"):
        vf = _open(kind, chrom)
        if vf is None:
            continue
        try:
            recs = list(vf.fetch(chrom, pos - 1, pos))
            # pysam parses record fields lazily, so accessing rec.info/.alts here
            # can raise a header-parse error ("Invalid header") on a flaky remote
            # read — keep the whole record loop inside the guard so it falls back
            # cleanly instead of propagating and crashing the caller.
            for rec in recs:
                if rec.pos != pos or rec.ref != variant.ref:
                    continue
                alts = rec.alts or ()
                if variant.alt not in alts:
                    continue
                cand = _best_from_record(rec)
                if cand and (best is None or (cand["af"] or 0) > (best["af"] or 0)):
                    best = cand
        except Exception:
            continue
        opened += 1
    if best is not None:
        return best            # a match from either callset is authoritative
    if opened == 2:
        # both callsets opened and neither carries the allele -> genuine absence.
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
    return None                # <2 callsets queryable -> can't assert absence, fall back
