"""PM1 (mutational hotspot) + PM5 at variable strength — both from the ClinVar residue index.

The residue index gives two distinct signals, and the split between them is what keeps ACMG's
"do not count PM1 with PM5" rule intact:

* **same residue** → PS1 (identical AA change) or PM5 (different AA change), graded by how well
  established the residue is;
* **neighbouring residues** → PM1, and only when the query's own residue carries nothing.
"""
import gzip

from vcf2report.acmg.criteria import _pm5_strength, all_criteria
from vcf2report.annotate import clinvar_residue
from vcf2report.models import Annotation, Variant

_pm1 = all_criteria()["PM1"]
_pm5 = all_criteria()["PM5"]


def _v(consequence="missense_variant", gene="TESTG", hgvs_p="p.Arg123Cys"):
    return Variant(chrom="1", pos=100, ref="C", alt="T", gene=gene,
                   consequence=consequence, hgvs_p=hgvs_p)


def _hot(n_residues=4, n_changes=9, window=7, enrichment=3.5):
    return {"n_residues": n_residues, "n_changes": n_changes, "window": window,
            "residues": [], "enrichment": enrichment, "gene_baseline": 0.08, "available": True}


def _ann(**kw):
    kw.setdefault("clinvar_residue_available", True)
    return Annotation(**kw)


# --- PM1 --------------------------------------------------------------------
def test_pm1_met_in_a_dense_neighbourhood():
    cr = _pm1(_v(), _ann(clinvar_hotspot=_hot()))
    assert cr.applies and cr.met
    assert cr.applied_strength == "moderate"
    assert cr.adjudicated_by == "engine"
    assert "hotspot" in cr.reasoning


def test_pm1_withheld_when_same_residue_evidence_exists():
    # ACMG: PM1 must not be counted together with PM5 — the same-residue evidence is PM5's.
    pm5_match = {"alt_aa": "His", "ref_aa": "Arg", "stars": 2, "n_other": 1}
    cr = _pm1(_v(), _ann(clinvar_hotspot=_hot(), clinvar_pm5=pm5_match))
    assert not cr.met
    assert "PS1/PM5" in cr.reasoning

    ps1_match = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2}
    assert not _pm1(_v(), _ann(clinvar_hotspot=_hot(), clinvar_ps1=ps1_match)).met


def test_pm1_not_met_below_the_density_cutoff():
    cr = _pm1(_v(), _ann(clinvar_hotspot=_hot(n_residues=2, n_changes=2)))
    assert not cr.met
    assert "needs" in cr.reasoning


def test_pm1_not_met_in_a_densely_catalogued_gene_without_local_enrichment():
    # The over-call guard: raw neighbour count alone fired on 61% of FBN1 / 49% of SCN1A novel
    # residues — exhaustively studied genes, not hotspots. A window must also be materially denser
    # than the gene's OWN baseline to count.
    cr = _pm1(_v(), _ann(clinvar_hotspot=_hot(n_residues=5, enrichment=1.1)))
    assert not cr.met
    assert "well-catalogued gene, not a local hotspot" in cr.reasoning
    assert cr.evidence["enrichment"] == 1.1


def test_pm1_not_met_in_a_missense_tolerant_gene():
    # Approximates ACMG's "without benign variation": a gene that tolerates missense is no hotspot.
    cr = _pm1(_v(), _ann(clinvar_hotspot=_hot(), gene_missense_tolerant=True))
    assert not cr.met
    assert "tolerates missense" in cr.reasoning


def test_pm1_not_met_for_non_missense():
    assert not _pm1(_v(consequence="stop_gained"), _ann(clinvar_hotspot=_hot())).met


def test_pm1_honest_when_index_unavailable():
    cr = _pm1(_v(), Annotation(clinvar_residue_available=False))
    assert not cr.met
    assert "unavailable" in cr.reasoning


# --- PM5 variable strength --------------------------------------------------
def test_pm5_strength_needs_both_breadth_and_review_quality():
    # Strong requires >=2 distinct changes AND a well-reviewed record. Breadth alone is not enough:
    # two single-submitter 1* records at one residue are suggestive, not Strong (that combination
    # otherwise carried a variant to Likely Pathogenic on two unreviewed submissions).
    assert _pm5_strength({"n_other": 3, "stars": 3}) == "strong"
    assert _pm5_strength({"n_other": 2, "stars": 1}) == "moderate"    # breadth, no quality
    assert _pm5_strength({"n_other": 1, "stars": 2}) == "moderate"    # quality, no breadth
    assert _pm5_strength({"n_other": 1, "stars": 1}) == "supporting"  # neither
    assert _pm5_strength(None) is None


def test_pm5_applies_the_graded_strength():
    strong = {"alt_aa": "His", "ref_aa": "Arg", "stars": 2, "n_other": 2,
              "accession": "VCV000001"}
    cr = _pm5(_v(), _ann(clinvar_pm5=strong))
    assert cr.met and cr.applied_strength == "strong"
    assert "PM5_strong" in cr.reasoning

    weak = {"alt_aa": "His", "ref_aa": "Arg", "stars": 1, "n_other": 1,
            "accession": "VCV000002"}
    cr = _pm5(_v(), _ann(clinvar_pm5=weak))
    assert cr.met and cr.applied_strength == "supporting"


def test_pp2_stands_down_when_pm1_fires():
    # PP2 and PM1 assert the same thing ("missense matters here") at gene vs. region granularity;
    # ClinGen VCEPs direct that they not both apply, so the more specific PM1 wins.
    _pp2 = all_criteria()["PP2"]
    a = _ann(clinvar_hotspot=_hot(), gene_mis_z=5.5, gene_missense_constrained=True)
    assert _pm1(_v(), a).met
    cr = _pp2(_v(), a)
    assert not cr.met
    assert "PM1" in cr.reasoning and "stands down" in cr.reasoning


def test_pp2_still_fires_without_a_hotspot():
    _pp2 = all_criteria()["PP2"]
    a = _ann(clinvar_hotspot=_hot(n_residues=0, enrichment=0.0),
             gene_mis_z=5.5, gene_missense_constrained=True)
    assert not _pm1(_v(), a).met
    assert _pp2(_v(), a).met


# --- hotspot() over a real-format table -------------------------------------
def _write_index(path, rows):
    with gzip.open(path, "wt") as w:
        w.write("# ClinVar residue index\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")


def test_hotspot_counts_neighbours_and_excludes_the_query_residue(tmp_path, monkeypatch):
    from vcf2report import config
    fp = tmp_path / "res.tsv.gz"
    _write_index(fp, [
        ("HOTG", 100, "Arg", "Cys", 2, "1-1-A-T", "VCV1"),   # the query residue itself
        ("HOTG", 102, "Gly", "Asp", 2, "1-2-A-T", "VCV2"),
        ("HOTG", 104, "Leu", "Pro", 1, "1-3-A-T", "VCV3"),
        ("HOTG", 105, "Val", "Met", 3, "1-4-A-T", "VCV4"),
        ("HOTG", 130, "Ser", "Asn", 2, "1-5-A-T", "VCV5"),   # outside the +/-7 window
    ])
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_FROZEN", fp)
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_LOCAL", None)
    clinvar_residue._index = None

    h = clinvar_residue.hotspot("HOTG", 100)
    assert h["n_residues"] == 3          # 102/104/105 — 100 itself excluded, 130 out of window
    assert h["residues"] == [102, 104, 105]
    assert h["n_changes"] == 3

    # a residue with no catalogued neighbours is not a hotspot
    assert clinvar_residue.hotspot("HOTG", 130)["n_residues"] == 0
    assert clinvar_residue.hotspot("NOSUCHGENE", 100)["n_residues"] == 0


def test_hotspot_enrichment_is_relative_to_the_genes_own_density(tmp_path, monkeypatch):
    """A tight cluster is enriched; the same raw count spread over a densely catalogued gene is not.

    Both genes below put 3 pathogenic residues inside the +-7 window of the query. CLUSTER has only
    those few catalogued in a long span (a genuine local concentration); DENSE has a pathogenic
    residue every 2 aa across the whole span (an exhaustively studied gene, no local signal).
    """
    from vcf2report import config
    rows = [("CLUSTER", p, "Arg", "Cys", 2, f"1-{p}-A-T", "VCV1") for p in (302, 304, 305)]
    rows += [("CLUSTER", 900, "Arg", "Cys", 2, "1-900-A-T", "VCV1")]      # stretches the span
    rows += [("DENSE", p, "Arg", "Cys", 2, f"2-{p}-A-T", "VCV2") for p in range(200, 400, 2)]
    fp = tmp_path / "res.tsv.gz"
    _write_index(fp, rows)
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_FROZEN", fp)
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_LOCAL", None)
    clinvar_residue._index = None

    tight = clinvar_residue.hotspot("CLUSTER", 303)
    spread = clinvar_residue.hotspot("DENSE", 303)
    assert tight["n_residues"] == 3 and spread["n_residues"] >= 3   # same raw density...
    assert tight["enrichment"] >= clinvar_residue.HOTSPOT_MIN_ENRICHMENT
    assert spread["enrichment"] < clinvar_residue.HOTSPOT_MIN_ENRICHMENT  # ...different meaning
