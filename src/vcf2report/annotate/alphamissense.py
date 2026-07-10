"""AlphaMissense missense-pathogenicity client (local tabix).

AlphaMissense (Cheng et al., Science 2023; (c) Google DeepMind, CC BY 4.0) scores
every possible human missense substitution 0..1. We read the tabix-indexed hg38
file (``AlphaMissense_hg38.tsv.gz``, fetched once by scripts/fetch_alphamissense.sh)
for a single variant — no network egress and no multi-GB data bundled in the repo.

File columns: ``#CHROM POS REF ALT genome uniprot_id transcript_id protein_variant
am_pathogenicity am_class``. A genomic variant can appear on several transcripts;
we take the maximum am_pathogenicity (most damaging) and its class. Non-missense
variants are simply absent from the file (score None) — LoF pathogenicity is PVS1's
job, not PP3's.
"""
from __future__ import annotations

import threading
from typing import Optional

from .. import config
from ..models import Variant
from . import cache

_SOURCE = "alphamissense"
_lock = threading.Lock()
_tabix = None
_tabix_tried = False
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


def _open():
    """Open (and cache) the local AlphaMissense TabixFile, or None if unavailable."""
    global _tabix, _tabix_tried
    if _tabix_tried:
        return _tabix
    with _lock:
        if _tabix_tried:
            return _tabix
        _tabix_tried = True
        pysam = _get_pysam()
        fp = config.ALPHAMISSENSE_LOCAL
        if pysam is None or not fp.exists():
            _tabix = None
            return None
        try:
            _tabix = pysam.TabixFile(str(fp))
        except Exception:
            _tabix = None
        return _tabix


def _fetch(tabix, chrom: str, pos: int) -> list[str]:
    """Rows overlapping (chrom, pos), trying with and without a 'chr' prefix."""
    bare = chrom[3:] if chrom.lower().startswith("chr") else chrom
    for c in (f"chr{bare}", bare):
        try:
            rows = list(tabix.fetch(c, pos - 1, pos))
        except Exception:
            continue
        if rows:
            return rows
    return []


def _best(rows: list[str], variant: Variant) -> Optional[dict]:
    """Max am_pathogenicity (and its class) among rows matching REF/ALT."""
    best_score = None
    best_class = None
    for row in rows:
        f = row.split("\t")
        if len(f) < 10:
            continue
        if f[2].upper() != variant.ref.upper() or f[3].upper() != variant.alt.upper():
            continue
        try:
            score = float(f[8])
        except (ValueError, IndexError):
            continue
        if best_score is None or score > best_score:
            best_score = score
            best_class = f[9].strip()
    if best_score is None:
        return None
    return {"am_pathogenicity": best_score, "am_class": best_class}


def lookup(variant: Variant) -> dict:
    """Return ``{'am_pathogenicity','am_class','_source'}`` for a variant.

    ``am_pathogenicity`` is None when AlphaMissense has no score (non-missense, or
    the local file / pysam is unavailable) — never a fabricated 0.
    """
    cached = cache.get(_SOURCE, variant.key)
    if cached is not None:
        return {**cached, "_source": "AlphaMissense (cache)"}

    tabix = _open()
    if tabix is None:
        return {"am_pathogenicity": None, "am_class": None,
                "_source": "AlphaMissense (unavailable — no local file/pysam)"}

    best = _best(_fetch(tabix, variant.chrom, variant.pos), variant)
    if best is None:
        result = {"am_pathogenicity": None, "am_class": None}
        cache.put(_SOURCE, variant.key, result)  # confirmed absence (non-missense)
        return {**result, "_source": "AlphaMissense hg38 (no missense score)"}
    cache.put(_SOURCE, variant.key, best)
    return {**best, "_source": "AlphaMissense hg38 (local tabix)"}
