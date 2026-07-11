"""Local gnomAD frequency client (reduced tabix table).

Reads a small, offline, tabix-indexed TSV built by ``scripts/build_gnomad_local.py``
that carries only the fields the ACMG engine cites — grpmax AF, AC/AN, homozygote
count, and the ClinGen filtering AF (faf95) — one row per gnomAD variant. It is the
same reduction ``gnomad_remote`` computes per record, so local and remote agree; the
local table just makes it offline and instant (no 150 GB full download, no network).

Absent file or missing pysam -> returns None so ``gnomad.lookup`` falls back to the
remote/live/bundled path exactly as before (behaviour-preserving).

Schema (tab-separated, bgzipped, ``tabix -s1 -b2 -e2``; ``#`` header skipped):
    #chrom  pos  ref  alt  af  ac  an  hom  faf95  pop
``chrom`` has no ``chr`` prefix (matches :pyattr:`Variant.key`); ``pop`` may be empty.
"""
from __future__ import annotations

import threading
from typing import Optional

from .. import config
from ..models import Variant

_tabix = None
_tabix_tried = False
_pysam = None
_pysam_tried = False
_mode: Optional[str] = None   # "full" | "partial"; drives the miss semantics
_lock = threading.Lock()


def _get_pysam():
    global _pysam, _pysam_tried
    if not _pysam_tried:
        _pysam_tried = True
        try:
            import pysam  # noqa: F401
            _pysam = pysam
        except Exception:
            _pysam = None
    return _pysam


def _open():
    """Lazily open the tabix file once; None if unavailable (missing file/pysam)."""
    global _tabix, _tabix_tried, _mode
    if _tabix_tried:
        return _tabix
    with _lock:
        if _tabix_tried:
            return _tabix
        _tabix_tried = True
        pysam = _get_pysam()
        fp = config.GNOMAD_LOCAL_TABIX
        if pysam is None or not fp.exists():
            _tabix = None
            return _tabix
        try:
            _tabix = pysam.TabixFile(str(fp))
        except Exception:
            _tabix = None
            return _tabix
        # Completeness: a full gnomAD build can assert absence on a miss; a partial
        # (panel / per-VCF) build cannot. Read it from a sidecar; default to the safe
        # "partial" so an unknown build never fabricates an absence.
        _mode = "partial"
        try:
            import json
            meta = fp.with_name(fp.name + ".meta")
            if meta.exists():
                m = (json.loads(meta.read_text()).get("mode") or "").lower()
                if m in ("full", "partial"):
                    _mode = m
        except Exception:
            pass
        return _tabix


def _norm_chrom(chrom: str) -> str:
    return chrom[3:] if chrom.lower().startswith("chr") else chrom


def _to_float(s: str) -> Optional[float]:
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def _to_int(s: str) -> Optional[int]:
    try:
        return int(s)
    except (TypeError, ValueError):
        return None


def query(variant: Variant) -> Optional[dict]:
    """Return {'af','ac','an','hom','faf95','pop'} for an exact CHROM-POS-REF-ALT
    match, ``{...af:0.0...}`` when the position is covered but this allele is absent,
    or None when no local table is available (caller falls back)."""
    tabix = _open()
    if tabix is None:
        return None
    chrom = _norm_chrom(variant.chrom)
    for name in (chrom, "chr" + chrom):   # tolerate either contig naming in the file
        try:
            rows = list(tabix.fetch(name, variant.pos - 1, variant.pos))
        except (ValueError, KeyError):
            rows = []
        if rows:
            break
    for row in rows:
        f = row.split("\t")
        if len(f) < 4:
            continue
        if _to_int(f[1]) != variant.pos or f[2] != variant.ref or f[3] != variant.alt:
            continue
        return {"af": _to_float(f[4]) if len(f) > 4 else None,
                "ac": _to_int(f[5]) if len(f) > 5 else None,
                "an": _to_int(f[6]) if len(f) > 6 else None,
                "hom": _to_int(f[7]) if len(f) > 7 else None,
                "faf95": _to_float(f[8]) if len(f) > 8 else None,
                "pop": (f[9].strip() or None) if len(f) > 9 else None}
    # No exact allele match. Only a FULL table (every gnomAD allele at every site)
    # may call that a true absence. A partial table (panel / per-VCF) is not even
    # position-complete — a --from-vcf table holds only the sample's own alleles — so
    # it returns None and the caller falls back, NEVER fabricating an absence.
    # (Genuinely-absent input variants still resolve: build writes them an explicit
    # af 0.0 row, which is matched above.)
    if _mode == "full":
        return {"af": 0.0, "ac": 0, "an": 0, "hom": 0, "faf95": 0.0, "pop": None}
    return None


def _reset_for_tests() -> None:
    global _tabix, _tabix_tried, _mode
    with _lock:
        _tabix, _tabix_tried, _mode = None, False, None
