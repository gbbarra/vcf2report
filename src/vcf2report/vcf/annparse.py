"""Parse consequence/HGVS from real annotator output in a VCF's INFO.

Supports SnpEff ``ANN`` and Ensembl VEP ``CSQ`` (whose subfield order is read
from the VCF header). Falls back to plain ``GENE``/``CSQ``/``HGVSC``/``HGVSP``
keys (as the bundled synthetic sample uses). Returns a dict with gene,
consequence, hgvs_c, hgvs_p, exon, transcript — or None if nothing usable is found.
"""
from __future__ import annotations

import re
from typing import Optional

# SnpEff ANN subfield order (fixed by the SnpEff spec).
_SNPEFF = ["allele", "annotation", "impact", "gene", "gene_id", "feature_type",
           "feature_id", "biotype", "rank", "hgvs_c", "hgvs_p"]


def parse_csq_format(header_lines: list[str]) -> Optional[list[str]]:
    """Extract the VEP CSQ subfield names from the ##INFO CSQ header."""
    for line in header_lines:
        if line.startswith("##INFO=<ID=CSQ") and "Format:" in line:
            fmt = line.split("Format:")[1].strip().rstrip('">').strip()
            return [f.strip() for f in fmt.split("|")]
    return None


def _first_term(consequence: str) -> Optional[str]:
    # VEP/SnpEff join multiple consequences with "&"; the first is most severe.
    return (consequence.split("&")[0].strip() or None) if consequence else None


def _minimal_alt(ref: str, alt: str) -> str:
    """VEP's minimal ALT: trim the longest shared leading base(s); '-' if emptied."""
    i = 0
    while i < len(ref) and i < len(alt) and ref[i] == alt[i]:
        i += 1
    return alt[i:] or "-"


def _allele_match(vep_allele: str, alt: str, ref: str = "") -> bool:
    if vep_allele == alt:
        return True
    if ref:
        return vep_allele == _minimal_alt(ref, alt)
    # No ref available: accept the single-base trim or a fully-trimmed insertion.
    return vep_allele == alt[1:] or vep_allele == "-"


def parse_snpeff(ann: str, alt: str, ref: str = "", n_alt: int = 1) -> Optional[dict]:
    first = None
    matched = None
    for entry in ann.split(","):
        f = entry.split("|")
        if len(f) < len(_SNPEFF):
            continue
        if first is None:
            first = f
        if _allele_match(f[0], alt, ref):
            matched = f
            break  # first allele-matching entry is most severe for that allele
    if matched is not None:
        f = matched
    elif n_alt == 1 and first is not None:
        f = first          # single-allele record: the only allele is safe to use
    else:
        return None        # multiallelic no-match: never borrow another allele
    return {"gene": f[3] or None, "consequence": _first_term(f[1]),
            "hgvs_c": f[9] or None, "hgvs_p": f[10] or None,
            "exon": f[8] or None,          # SnpEff "rank" = exon "N/M"
            "transcript": f[6] or None}    # SnpEff "feature_id" = transcript


def parse_vep(csq: str, alt: str, field_names: list[str], ref: str = "",
              alt_index: int = 0, n_alt: int = 1) -> Optional[dict]:
    from urllib.parse import unquote
    idx = {name.lower(): i for i, name in enumerate(field_names)}

    def get(f: list[str], name: str) -> Optional[str]:
        i = idx.get(name)
        v = f[i] if (i is not None and i < len(f) and f[i]) else None
        return unquote(v) if v else v  # VEP percent-encodes reserved chars

    def row(f: list[str]) -> dict:
        return {"gene": get(f, "symbol") or get(f, "gene"),
                "consequence": _first_term(get(f, "consequence") or ""),
                "hgvs_c": get(f, "hgvsc"), "hgvs_p": get(f, "hgvsp"),
                "exon": get(f, "exon"),
                "transcript": get(f, "feature")}   # VEP "Feature" = transcript id

    entries = [e.split("|") for e in csq.split(",")]
    an = idx.get("allele_num")
    a = idx.get("allele")

    def matches(f: list[str]) -> bool:
        if an is not None and an < len(f) and f[an]:   # prefer VEP's ALLELE_NUM
            return f[an] == str(alt_index + 1)
        return a is not None and a < len(f) and _allele_match(f[a], alt, ref)

    candidates = [f for f in entries if matches(f)]
    if not candidates:
        if n_alt == 1:
            candidates = entries          # single allele: all blocks are this ALT
        else:
            return None                   # multiallelic no-match: don't borrow
    for f in candidates:                  # PICK > CANONICAL > MANE_SELECT > first
        if get(f, "pick") == "1":
            return row(f)
    for f in candidates:
        if get(f, "canonical") == "YES":
            return row(f)
    for f in candidates:
        if get(f, "mane_select"):
            return row(f)
    return row(candidates[0])


# Legacy SnpEff EFF effect names -> Sequence Ontology terms the engine uses.
_EFF_SO = {
    "NON_SYNONYMOUS_CODING": "missense_variant", "SYNONYMOUS_CODING": "synonymous_variant",
    "STOP_GAINED": "stop_gained", "STOP_LOST": "stop_lost", "START_LOST": "start_lost",
    "NON_SYNONYMOUS_START": "start_lost", "FRAME_SHIFT": "frameshift_variant",
    "SPLICE_SITE_ACCEPTOR": "splice_acceptor_variant", "SPLICE_SITE_DONOR": "splice_donor_variant",
    "CODON_DELETION": "inframe_deletion", "CODON_INSERTION": "inframe_insertion",
    "CODON_CHANGE_PLUS_CODON_DELETION": "inframe_deletion",
    "CODON_CHANGE_PLUS_CODON_INSERTION": "inframe_insertion", "EXON_DELETED": "transcript_ablation",
    "SYNONYMOUS_STOP": "stop_retained_variant", "UTR_5_PRIME": "5_prime_UTR_variant",
    "UTR_3_PRIME": "3_prime_UTR_variant", "INTRON": "intron_variant",
    "UPSTREAM": "upstream_gene_variant", "DOWNSTREAM": "downstream_gene_variant",
    "INTERGENIC": "intergenic_variant",
}
_EFF_IMPACT = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
_AA3 = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
        "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
        "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "*": "Ter"}


def _eff_hgvs_p(aa: str) -> Optional[str]:
    """SnpEff EFF one-letter AA change ('W169R', 'Q59*') -> HGVS.p ('p.Trp169Arg')."""
    m = re.match(r"^([A-Z*])(\d+)([A-Z*])$", aa or "")
    if not m:
        return None
    a, pos, b = m.groups()
    return f"p.{_AA3.get(a, a)}{pos}{_AA3.get(b, b)}"


def parse_snpeff_eff(eff: str, alt_index: int = 0, n_alt: int = 1) -> Optional[dict]:
    """Parse the legacy SnpEff ``EFF`` INFO (pre-``ANN`` format, still common in older
    public VCFs). Picks the most severe effect (impact, protein-coding preferred) for
    this ALT; EFF carries no HGVS.c so that stays None."""
    best, best_score = None, -1
    for entry in eff.split(","):
        if "(" not in entry:
            continue
        effect, rest = entry.split("(", 1)
        f = rest.rstrip(")").split("|")
        if len(f) < 10:
            continue
        gt = f[10] if len(f) > 10 else ""
        if n_alt > 1 and gt.isdigit() and int(gt) != alt_index + 1:
            continue  # multiallelic: only this ALT's effects
        score = _EFF_IMPACT.get(f[0], 0) * 2 + (1 if f[6] == "protein_coding" else 0)
        if score > best_score:
            best_score, best = score, (effect.strip(), f)
    if not best:
        return None
    effect, f = best
    return {"gene": f[5] or None,
            "consequence": _EFF_SO.get(effect, effect.lower() or None),
            "hgvs_c": None, "hgvs_p": _eff_hgvs_p(f[3]),
            "exon": f[9] or None, "transcript": f[8] or None}


def extract(info: dict[str, str], alt: str, csq_format: Optional[list[str]] = None,
            ref: str = "", alt_index: int = 0, n_alt: int = 1) -> Optional[dict]:
    """Best consequence/HGVS from ANN, then legacy EFF, then CSQ, then plain keys."""
    if info.get("ANN"):
        r = parse_snpeff(info["ANN"], alt, ref, n_alt)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    if info.get("EFF"):
        r = parse_snpeff_eff(info["EFF"], alt_index, n_alt)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    csq = info.get("CSQ")
    # Only treat CSQ as VEP if it's the structured, pipe-delimited form.
    if csq and "|" in csq and csq_format:
        r = parse_vep(csq, alt, csq_format, ref, alt_index, n_alt)
        if r and (r.get("gene") or r.get("consequence")):
            return r
    # Plain keys (synthetic sample / simple pipelines).
    simple = {"gene": info.get("GENE"),
              "consequence": info.get("CSQ") if (info.get("CSQ") and "|" not in info["CSQ"]) else None,
              "hgvs_c": info.get("HGVSC"), "hgvs_p": info.get("HGVSP"),
              "transcript": info.get("TRANSCRIPT") or info.get("FEATURE")}
    return simple if any(simple.values()) else None
