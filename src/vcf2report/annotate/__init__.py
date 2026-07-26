"""Annotation aggregation: merge every source into one :class:`Annotation`.

``annotate_variant`` is the single entry point used by the headless pipeline and
the MCP ``annotate`` tool. It records, per field, which source (and date) it came
from so the ACMG criteria and the report can cite provenance.
"""
from __future__ import annotations

from ..models import Annotation, Variant
from . import abraom, alphamissense, clinvar, clinvar_residue, extra, from_vcf, gnomad, hpo


def annotate_variant(variant: Variant, patient_hpo: list[str] | None = None,
                     build_trusted: bool = True,
                     with_alphamissense: bool = True,
                     with_clinvar_residue: bool = True) -> Annotation:
    """Merge all sources into an :class:`Annotation`.

    ``build_trusted=False`` (a detected genome-build mismatch) means the variant's
    coordinates cannot be trusted against the GRCh38 databases, so all
    coordinate-keyed lookups are skipped and gnomAD AF is left unknown (None) —
    this prevents PM2 firing on a cross-build "absent" that is really a
    wrong-position miss. Gene-level (constraint, HPO) evidence stays valid.

    ``with_alphamissense=False`` skips only the AlphaMissense *client* lookup (a
    per-variant tabix hit on a ~1 GB file) — a pre-annotated INFO score is still
    read. The pipeline uses this to annotate the whole post-QC set cheaply, then
    enriches just the surviving candidates via :func:`add_alphamissense` (the score
    only feeds PP3/BP4 at classification time, never the filter, so deferring it is
    behaviour-preserving for the classified variants).

    ``with_clinvar_residue=False`` likewise skips the residue-index lookup that feeds
    PS1/PM5/PM1, for the pipeline to defer to :func:`add_clinvar_residue` on the
    candidates. It defaults to TRUE so that every ad-hoc caller — the MCP
    ``classify_variant`` tool, scripts, notebooks — gets a COMPLETE annotation from this
    one call. Making the deferral unconditional silently stripped those three criteria
    from the MCP path, which then told users to build an index they already had.
    """
    patient_hpo = patient_hpo or []

    if build_trusted:
        # Prefer annotations already in the VCF INFO (SnpEff/VEP + vcfanno) — the
        # fast, offline path for a real pre-annotated exome — then fall back to the
        # local snapshots / live clients only for whatever INFO didn't provide.
        vi = from_vcf.extract(variant)
        if "gnomad_af" in vi:
            g = {"af": vi["gnomad_af"], "ac": vi.get("gnomad_ac"),
                 "an": vi.get("gnomad_an"), "hom": vi.get("gnomad_hom"),
                 "faf95": vi.get("gnomad_faf95"), "pop": None, "_source": "VCF INFO"}
        else:
            g = gnomad.lookup(variant)
        if "clinvar_significance" in vi:
            cv = {"significance": vi["clinvar_significance"],
                  "review_status": vi.get("clinvar_review_status"),
                  "accession": vi.get("clinvar_accession"),
                  "condition": vi.get("clinvar_condition"), "date": None,
                  "_source": "VCF INFO"}
        else:
            cv = clinvar.lookup(variant)
        ab = {"af": vi["abraom_af"], "_source": "VCF INFO"} if "abraom_af" in vi \
            else abraom.lookup(variant)
        if "revel" in vi or "cadd" in vi:
            isi = {"revel": vi.get("revel"), "cadd": vi.get("cadd"), "_source": "VCF INFO"}
        else:
            isi = extra.insilico(variant)
        if "am_pathogenicity" in vi:
            am = {"am_pathogenicity": vi["am_pathogenicity"],
                  "am_class": vi.get("am_class"), "_source": "VCF INFO"}
        elif with_alphamissense:
            am = alphamissense.lookup(variant)
        else:
            am = {"am_pathogenicity": None, "am_class": None,
                  "_source": "AlphaMissense (deferred to candidate stage)"}
    else:
        note = "skipped — genome-build mismatch (coordinates not GRCh38)"
        g = {"af": None, "ac": None, "an": None, "hom": None, "pop": None, "_source": note}
        cv = {"significance": None, "review_status": None, "accession": None,
              "condition": None, "date": None, "_source": note}
        ab = {"af": None, "_source": note}
        isi = {"revel": None, "cadd": None, "_source": note}
        am = {"am_pathogenicity": None, "am_class": None, "_source": note}
    con = extra.gene_constraint(variant.gene)
    ph = hpo.match(variant.gene, patient_hpo)

    ann = Annotation(
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
        gnomad_faf95=g.get("faf95"),
        abraom_af=ab.get("af"),
        gene_lof_intolerant=con.get("lof_intolerant"),
        gene_mis_z=con.get("mis_z"),
        gene_oe_mis_upper=con.get("oe_mis_upper"),
        gene_missense_constrained=con.get("missense_constrained"),
        gene_missense_tolerant=con.get("missense_tolerant"),
        revel=isi.get("revel"),
        cadd_phred=isi.get("cadd"),
        am_pathogenicity=am.get("am_pathogenicity"),
        am_class=am.get("am_class"),
        hpo_match_score=ph.get("score"),
        hpo_best_match=ph.get("best"),
        hpo_matched_terms=ph.get("matched_terms", []),
        source={
            "gnomad": g.get("_source", ""),
            "clinvar": cv.get("_source", ""),
            "abraom": ab.get("_source", ""),
            # One provenance string for the whole gnomAD constraint table — PVS1 cites it for the
            # LoF columns, PP2/BP1 for the missense ones. A single key, so the two can never drift.
            "gene_constraint": con.get("_source", ""),
            "insilico": isi.get("_source", ""),
            "alphamissense": am.get("_source", ""),
            "hpo": ph.get("_source", ""),
        },
    )
    if with_clinvar_residue:
        add_clinvar_residue(variant, ann, build_trusted=build_trusted)
    return ann


def add_clinvar_residue(variant: Variant, annotation: Annotation,
                        build_trusted: bool = True) -> None:
    """Populate the residue-level ClinVar evidence (PS1/PM5/PM1) — the lazy path.

    Deferred for the same reason as :func:`add_alphamissense`: these fields feed classification
    ONLY (``ps1``/``pm5``/``pm1``/``pp2``), never the QC/rarity filter, so computing them for every
    post-QC variant does ~10-20k lookups per exome to serve the few hundred that survive to
    candidacy. The pipeline calls this on the candidates instead, which is behaviour-preserving.

    Needs the GRCh38 genomic key to prove a PS1 hit is a *different* variant, so it is gated on a
    trusted build like the other coordinate-keyed lookups; a build mismatch leaves the matches
    unpopulated (None) rather than fabricating a cross-build match.
    """
    if not build_trusted:
        return
    res = clinvar_residue.lookup(variant.gene, variant.hgvs_p, variant.key)
    annotation.clinvar_ps1 = res.get("ps1")
    annotation.clinvar_pm5 = res.get("pm5")
    annotation.clinvar_residue_available = res.get("available")
    # `lookup` already parsed the protein position — reuse it rather than re-running the regex.
    aa_pos = res.get("aa_pos")
    if aa_pos is not None:
        annotation.clinvar_hotspot = clinvar_residue.hotspot(variant.gene, aa_pos)
    annotation.source["clinvar_residue"] = "ClinVar residue index (local)"


def add_alphamissense(variant: Variant, annotation: Annotation) -> None:
    """Populate AlphaMissense on an already-built Annotation (the lazy path).

    Called for the surviving candidates only, so the ~1 GB tabix file is queried a
    few times per report instead of once per post-QC variant. A no-op if the score
    is already present (e.g. read from VCF INFO). The lookup is identical to the one
    ``annotate_variant`` would have run, so the resulting classification is unchanged.
    """
    if annotation.am_pathogenicity is not None:
        return
    am = alphamissense.lookup(variant)
    annotation.am_pathogenicity = am.get("am_pathogenicity")
    annotation.am_class = am.get("am_class")
    annotation.source["alphamissense"] = am.get("_source", "")
