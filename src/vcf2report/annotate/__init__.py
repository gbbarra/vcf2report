"""Annotation aggregation: merge every source into one :class:`Annotation`.

``annotate_variant`` is the single entry point used by the headless pipeline and
the MCP ``annotate`` tool. It records, per field, which source (and date) it came
from so the ACMG criteria and the report can cite provenance.
"""
from __future__ import annotations

from ..models import Annotation, Variant
from . import abraom, clinvar, extra, gnomad, hpo


def annotate_variant(variant: Variant, patient_hpo: list[str] | None = None) -> Annotation:
    patient_hpo = patient_hpo or []

    g = gnomad.lookup(variant)
    cv = clinvar.lookup(variant)
    ab = abraom.lookup(variant)
    con = extra.gene_constraint(variant.gene)
    isi = extra.insilico(variant)
    ph = hpo.match(variant.gene, patient_hpo)

    return Annotation(
        clinvar_significance=cv.get("significance"),
        clinvar_review_status=cv.get("review_status"),
        clinvar_accession=cv.get("accession"),
        clinvar_condition=cv.get("condition"),
        clinvar_date=cv.get("date"),
        gnomad_af=g.get("af"),
        gnomad_ac=g.get("ac"),
        gnomad_an=g.get("an"),
        gnomad_homozygotes=g.get("hom"),
        gnomad_popmax_pop=g.get("pop"),
        abraom_af=ab.get("af"),
        gene_lof_intolerant=con.get("lof_intolerant"),
        revel=isi.get("revel"),
        cadd_phred=isi.get("cadd"),
        hpo_match_score=ph.get("score"),
        hpo_matched_terms=ph.get("matched_terms", []),
        source={
            "gnomad": g.get("_source", ""),
            "clinvar": cv.get("_source", ""),
            "abraom": ab.get("_source", ""),
            "gene_lof_intolerant": con.get("_source", ""),
            "insilico": isi.get("_source", ""),
            "hpo": ph.get("_source", ""),
        },
    )
