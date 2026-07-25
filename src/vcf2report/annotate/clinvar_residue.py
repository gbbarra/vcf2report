"""ClinVar residue index — the lookup behind PS1 and PM5.

PS1 (same amino-acid change as a *different* established pathogenic variant) and PM5
(novel missense at a residue where a *different* missense is pathogenic) both need to
know, for a query missense ``(gene, protein_position, alt_aa)``, which pathogenic
amino-acid changes ClinVar has already reported at that residue.

The table (built by ``scripts/fetch_clinvar_residue.py``) has one row per distinct
``(gene, aa_pos, alt_aa)`` pathogenic/likely-pathogenic missense with a criteria-based
(>=1-star) review:

    gene<TAB>aa_pos<TAB>ref_aa<TAB>alt_aa<TAB>stars<TAB>genomic_key<TAB>accession

Two sources are merged, both optional: the committed *frozen slice*
(``CLINVAR_RESIDUE_FROZEN``, a small demo/test subset) and the full genome-wide index
(``CLINVAR_RESIDUE_LOCAL``, built locally, git-ignored). The local index wins on
conflicts. When neither is present the lookup returns empty matches, so PS1/PM5 report
"index unavailable" rather than a fabricated hit.
"""
from __future__ import annotations

import re
from typing import Optional

from .. import config

# gene -> {aa_pos(int) -> {alt_aa -> (ref_aa, stars, genomic_key, accession)}}
_index: Optional[dict] = None

_P_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
_AA3 = {"Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
        "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"}


def parse_hgvs_p(hgvs_p: Optional[str]) -> Optional[tuple[str, int, str]]:
    """``p.Ser330Asn`` -> ``("Ser", 330, "Asn")`` for a clean missense, else None.

    Only substitutions between two standard amino acids qualify; stop/synonymous/indel
    ``p.`` forms (Ter, ``=``, dup, del, fs) return None so PS1/PM5 never fire on them.
    """
    if not hgvs_p:
        return None
    m = _P_RE.search(hgvs_p)
    if not m:
        return None
    ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3)
    if ref_aa not in _AA3 or alt_aa not in _AA3 or ref_aa == alt_aa:
        return None
    return ref_aa, pos, alt_aa


def _read_rows(fp):
    if not fp or not fp.exists():
        return
    import gzip
    with gzip.open(fp, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("gene\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            gene, pos, ref_aa, alt_aa, stars = p[0], p[1], p[2], p[3], p[4]
            key = p[5] if len(p) > 5 else None
            acc = p[6] if len(p) > 6 else None
            try:
                yield gene, int(pos), ref_aa, alt_aa, int(stars), key, acc
            except ValueError:
                continue


def _load() -> dict:
    global _index
    if _index is None:
        d: dict = {}
        # frozen slice first, then the local full index (which overrides on conflict).
        for fp in (config.CLINVAR_RESIDUE_FROZEN, config.CLINVAR_RESIDUE_LOCAL):
            for gene, pos, ref_aa, alt_aa, stars, key, acc in _read_rows(fp):
                residues = d.setdefault(gene, {}).setdefault(pos, {})
                prev = residues.get(alt_aa)
                if prev is None or stars >= prev[1]:
                    residues[alt_aa] = (ref_aa, stars, key, acc)
        _index = d  # publish only when fully built
    return _index


def available() -> bool:
    return bool(_load())


# PM1 hotspot window: distinct pathogenic-missense RESIDUES within +/- this many residues of the
# query (the query's own residue excluded — that is PS1/PM5 territory, never PM1).
HOTSPOT_WINDOW = 7
HOTSPOT_MIN_RESIDUES = 3
# ...AND the window must be this many times denser than the gene's OWN baseline density.
# Raw density alone measures how thoroughly a gene has been STUDIED, not where its hot spots are:
# on the committed slice a bare >=3-residue rule fired on 61% of FBN1 and 49% of SCN1A — exhaustively
# catalogued genes where almost every position has three pathogenic neighbours. Requiring local
# enrichment over the gene's own average restores the actual meaning of "hot spot": a local
# CONCENTRATION, not a well-covered gene.
HOTSPOT_MIN_ENRICHMENT = 2.0


def hotspot(gene: Optional[str], aa_pos: Optional[int], window: int = HOTSPOT_WINDOW) -> dict:
    """Local density of pathogenic missense around a residue — the PM1 signal.

    Counts DISTINCT neighbouring residues (excluding ``aa_pos`` itself) that carry at least one
    P/LP missense in the ClinVar residue index, plus the total pathogenic amino-acid changes across
    them, and expresses that as an ``enrichment`` over the gene's own baseline density (catalogued
    pathogenic residues per residue of its catalogued span). The query residue is excluded so PM1
    never double-counts the same-residue evidence PS1 and PM5 already carry.
    """
    idx = _load()
    out = {"n_residues": 0, "n_changes": 0, "window": window, "residues": [],
           "enrichment": 0.0, "gene_baseline": 0.0, "available": bool(idx)}
    if not gene or aa_pos is None:
        return out
    residues = idx.get(gene)
    if not residues:
        return out
    near = [p for p in residues
            if p != aa_pos and abs(p - aa_pos) <= window and residues[p]]
    out["n_residues"] = len(near)
    out["n_changes"] = sum(len(residues[p]) for p in near)
    out["residues"] = sorted(near)

    # gene baseline: catalogued pathogenic residues per residue across the catalogued span.
    lo, hi = min(residues), max(residues)
    span = max(hi - lo + 1, 1)
    baseline = len(residues) / span
    local = len(near) / float(2 * window + 1)
    out["gene_baseline"] = round(baseline, 4)
    out["enrichment"] = round(local / baseline, 2) if baseline else 0.0
    return out


def lookup(gene: Optional[str], hgvs_p: Optional[str], variant_key: Optional[str]) -> dict:
    """Residue matches for PS1 / PM5.

    Returns ``{"ps1": match|None, "pm5": match|None, "available": bool, "residue": str|None}``
    where a match is ``{"alt_aa","ref_aa","stars","accession","genomic_key"}``.

    * **ps1**: a ClinVar P/LP missense with the *same* amino-acid change at a *different*
      genomic locus (a distinct variant — the query's own record is PP5, not PS1).
    * **pm5**: a ClinVar P/LP missense with a *different* amino-acid change at the same
      residue, applied only when the query's exact change is *not itself* established
      (so PS1 and PM5 are mutually exclusive).
    """
    idx = _load()
    out = {"ps1": None, "pm5": None, "available": bool(idx), "residue": None}
    parsed = parse_hgvs_p(hgvs_p)
    if not gene or parsed is None:
        return out
    ref_aa, pos, alt_aa = parsed
    out["residue"] = f"{ref_aa}{pos}"
    residues = idx.get(gene, {}).get(pos)
    if not residues:
        return out

    # PS1: same amino-acid change, established pathogenic, at a DIFFERENT variant.
    same = residues.get(alt_aa)
    known_same = same is not None
    if known_same and same[2] and variant_key and same[2] != variant_key:
        out["ps1"] = {"alt_aa": alt_aa, "ref_aa": same[0], "stars": same[1],
                      "genomic_key": same[2], "accession": same[3]}
    elif known_same and (same[2] is None or not variant_key):
        # Same AA change is established but we can't prove it is a *different* variant
        # (missing key). Treat conservatively as PS1 evidence — a same-AA pathogenic exists.
        out["ps1"] = {"alt_aa": alt_aa, "ref_aa": same[0], "stars": same[1],
                      "genomic_key": same[2], "accession": same[3]}

    # PM5: a DIFFERENT pathogenic AA change at this residue — only when the query's exact
    # change is not itself established (else it is PS1/PP5 territory, never PM5). ``n_other``
    # counts the distinct pathogenic changes at the residue, which grades PM5's strength.
    if not known_same:
        best, n_other = None, 0
        for other_alt, (o_ref, o_stars, o_key, o_acc) in residues.items():
            if other_alt == alt_aa:
                continue
            n_other += 1
            if best is None or o_stars > best["stars"]:
                best = {"alt_aa": other_alt, "ref_aa": o_ref, "stars": o_stars,
                        "genomic_key": o_key, "accession": o_acc}
        if best is not None:
            best["n_other"] = n_other
        out["pm5"] = best
    return out
