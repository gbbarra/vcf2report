"""Read a GA4GH Phenopacket (v2 JSON) into pipeline inputs.

Phenopackets are the literature-endorsed standard for carrying phenotype + genotype
together (GA4GH). Public repositories (e.g. the GA4GH/Monarch Phenopacket Store)
hold thousands of curated, published patient cases — a practical source of paired
HPO terms + causative variants for testing vcf2report end to end.

This extracts the HPO term ids (excluded features skipped) and the variant(s)
(coordinates from ``vcfRecord``, gene from ``geneContext``, HGVS from
``expressions``, zygosity from ``allelicState``). Molecular *consequence* is not
part of a phenopacket, so the emitted VCF should be run through the annotation
recipe (docs/ANNOTATION.md) before classification — unless it is already present.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

_ZYGOSITY = {
    "GENO:0000135": "het", "heterozygous": "het",
    "GENO:0000136": "hom", "homozygous": "hom",
    "GENO:0000134": "hemi", "hemizygous": "hemi",
}


def _hpo_terms(pkt: dict) -> list[str]:
    terms: list[str] = []
    for feat in pkt.get("phenotypicFeatures", []) or []:
        if feat.get("excluded"):
            continue
        tid = (feat.get("type") or {}).get("id")
        if tid and tid.startswith("HP:"):
            terms.append(tid)
    return terms


def _zygosity_of(descriptor: dict) -> str | None:
    state = descriptor.get("allelicState") or {}
    return _ZYGOSITY.get(state.get("id")) or _ZYGOSITY.get((state.get("label") or "").lower())


def _variant_from_descriptor(descriptor: dict) -> dict | None:
    rec = descriptor.get("vcfRecord") or {}
    # Require full coordinates + alleles; an HGVS-only descriptor cannot be emitted
    # as a VCF row without resolving HGVS to coordinates.
    if not rec.get("chrom") or rec.get("pos") in (None, "") \
            or not rec.get("ref") or not rec.get("alt"):
        return None
    exprs = {e.get("syntax"): e.get("value") for e in descriptor.get("expressions", []) or []}
    gene = (descriptor.get("geneContext") or {}).get("symbol")
    return {
        "chrom": str(rec["chrom"]),
        "pos": int(rec["pos"]),
        "ref": rec.get("ref", ""),
        "alt": rec.get("alt", ""),
        "gene": gene,
        "hgvs_c": exprs.get("hgvs.c"),
        "hgvs_p": exprs.get("hgvs.p"),
        "zygosity": _zygosity_of(descriptor),
    }


def _variants(pkt: dict) -> tuple[list[dict], int]:
    """Return (emittable variants, count of genotype entries skipped)."""
    out: list[dict] = []
    skipped = 0
    for interp in pkt.get("interpretations", []) or []:
        diag = interp.get("diagnosis") or {}
        for gi in diag.get("genomicInterpretations", []) or []:
            vi = gi.get("variantInterpretation") or {}
            desc = vi.get("variationDescriptor") or {}
            v = _variant_from_descriptor(desc)
            if v:
                out.append(v)
            elif desc:  # a descriptor was present but lacked usable coordinates
                skipped += 1
    return out, skipped


def _disease(pkt: dict) -> str:
    """The case diagnosis label. Prefer an interpretation's ``diagnosis.disease``;
    fall back to a top-level ``diseases[]`` term. Empty string if none is recorded.
    (Phenopackets carry the diagnosis here — it is NOT part of ``vcfRecord`` — so a
    cohort builder must read it from the packet, not re-derive it from the variant.)"""
    for interp in pkt.get("interpretations", []) or []:
        label = ((interp.get("diagnosis") or {}).get("disease") or {}).get("label")
        if label:
            return label
    for dz in pkt.get("diseases", []) or []:
        if dz.get("excluded"):
            continue
        label = (dz.get("term") or {}).get("label")
        if label:
            return label
    return ""


def load_phenopacket(path: str | Path) -> dict[str, Any]:
    """Return {'subject_id', 'hpo_terms', 'disease', 'variants', 'skipped_variants'}."""
    pkt = json.loads(Path(path).read_text())
    variants, skipped = _variants(pkt)
    return {
        "subject_id": (pkt.get("subject") or {}).get("id") or "PHENOPACKET",
        "hpo_terms": _hpo_terms(pkt),
        "disease": _disease(pkt),
        "variants": variants,
        "skipped_variants": skipped,
    }


_VCF_HEADER = """##fileformat=VCFv4.2
##reference=GRCh38
##source=vcf2report_phenopacket
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=HGVSC,Number=1,Type=String,Description="HGVS coding">
##INFO=<ID=HGVSP,Number=1,Type=String,Description="HGVS protein">
##comment=From a GA4GH Phenopacket. Molecular consequence is NOT included — annotate
##comment=with docs/ANNOTATION.md (SnpEff) before classification unless already present.
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

_GT = {"het": "0/1", "hom": "1/1", "hemi": "1"}


def write_inputs(data: dict, vcf_path: str | Path, hpo_path: str | Path) -> None:
    """Write a VCF + an HPO-terms file from load_phenopacket() output."""
    vcf_path, hpo_path = Path(vcf_path), Path(hpo_path)
    def esc(x: str) -> str:
        # INFO uses ';' (field sep), '=' (key/val), and whitespace is illegal.
        return str(x).replace(";", "|").replace("=", "_").replace(" ", "_").replace("\t", "_")

    lines = [_VCF_HEADER.format(sample=data["subject_id"]).rstrip("\n")]
    for v in data["variants"]:
        chrom = v["chrom"][3:] if v["chrom"].lower().startswith("chr") else v["chrom"]
        info = ";".join(p for p in [
            f"GENE={esc(v['gene'])}" if v.get("gene") else "",
            f"HGVSC={esc(v['hgvs_c'])}" if v.get("hgvs_c") else "",
            f"HGVSP={esc(v['hgvs_p'])}" if v.get("hgvs_p") else "",
        ] if p) or "."
        gt = _GT.get(v.get("zygosity") or "het", "0/1")
        lines.append("\t".join([chrom, str(v["pos"]), ".", v["ref"], v["alt"],
                                "800", "PASS", info, "GT", gt]))
    vcf_path.write_text("\n".join(lines) + "\n")
    hpo_path.write_text("\n".join(data["hpo_terms"]) + "\n")
