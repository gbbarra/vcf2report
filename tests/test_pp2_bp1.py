"""PP2 / BP1 — gene-level missense-constraint criteria (now engine-decided).

Both read gnomAD per-gene missense constraint baked into the Annotation:

* **PP2** (Supporting pathogenic): a missense variant in a gene significantly
  *depleted* of missense variation (``mis_z`` >= 3.09).
* **BP1** (Supporting benign): a missense variant in a gene that is LoF-intolerant
  (truncating is the disease mechanism) yet *tolerates* missense
  (``oe_mis_upper`` >= 1.0).

They are mutually exclusive by construction — a gene cannot be both depleted of and
tolerant to missense — so no variant ever earns PP2 and BP1 together.
"""
from vcf2report.acmg.criteria import all_criteria
from vcf2report.annotate import extra
from vcf2report.models import Annotation, Variant

_pp2 = all_criteria()["PP2"]
_bp1 = all_criteria()["BP1"]


def _v(consequence="missense_variant", gene="TESTG"):
    return Variant(chrom="1", pos=100, ref="A", alt="G", gene=gene,
                   consequence=consequence, zygosity="het")


# --- PP2 --------------------------------------------------------------------
def test_pp2_met_missense_in_constrained_gene():
    a = Annotation(gene_mis_z=5.5, gene_missense_constrained=True,
                   source={"gene_constraint": "gnomAD v2.1.1 constraint (local)"})
    cr = _pp2(_v(), a)
    assert cr.applies and cr.met
    assert cr.applied_strength == "supporting"
    assert cr.adjudicated_by == "engine"
    assert cr.citation == ["gnomAD v2.1.1 constraint (local)"]


def test_pp2_not_met_for_non_missense():
    a = Annotation(gene_mis_z=5.5, gene_missense_constrained=True)
    assert not _pp2(_v(consequence="stop_gained"), a).met


def test_pp2_not_met_when_gene_not_constrained():
    a = Annotation(gene_mis_z=1.2, gene_missense_constrained=False)
    cr = _pp2(_v(), a)
    assert cr.applies and not cr.met
    assert "below" in cr.reasoning


def test_pp2_not_met_when_no_metric():
    cr = _pp2(_v(), Annotation())
    assert cr.applies and not cr.met
    assert "no gnomAD missense-constraint metric" in cr.reasoning


# --- BP1 --------------------------------------------------------------------
def test_bp1_met_missense_in_lof_intolerant_tolerant_gene():
    a = Annotation(gene_lof_intolerant=True, gene_oe_mis_upper=1.08,
                   gene_missense_tolerant=True, gene_missense_constrained=False,
                   source={"gene_constraint": "gnomAD v2.1.1 constraint (local)"})
    cr = _bp1(_v(), a)
    assert cr.applies and cr.met
    assert cr.applied_strength == "supporting"
    assert cr.adjudicated_by == "engine"


def test_bp1_not_met_for_non_missense():
    a = Annotation(gene_lof_intolerant=True, gene_missense_tolerant=True)
    assert not _bp1(_v(consequence="stop_gained"), a).met


def test_bp1_not_met_when_gene_not_lof_intolerant():
    # Missense-tolerant but not LoF-intolerant → truncating-only mechanism not established.
    a = Annotation(gene_lof_intolerant=False, gene_missense_tolerant=True)
    cr = _bp1(_v(), a)
    assert not cr.met
    assert "not LoF-intolerant" in cr.reasoning


def test_bp1_not_met_when_gene_missense_constrained():
    # A missense-constrained gene is never BP1, even if oe_mis flag were somehow set.
    a = Annotation(gene_lof_intolerant=True, gene_missense_tolerant=True,
                   gene_missense_constrained=True)
    assert not _bp1(_v(), a).met


def test_pp2_bp1_mutually_exclusive():
    # PP2's constrained gate and BP1's tolerant gate can never both be True.
    constrained = Annotation(gene_mis_z=5.5, gene_missense_constrained=True,
                             gene_lof_intolerant=True, gene_missense_tolerant=False)
    tolerant = Annotation(gene_lof_intolerant=True, gene_missense_tolerant=True,
                          gene_missense_constrained=False, gene_oe_mis_upper=1.08)
    assert _pp2(_v(), constrained).met and not _bp1(_v(), constrained).met
    assert _bp1(_v(), tolerant).met and not _pp2(_v(), tolerant).met


# --- integration through the real local constraint table --------------------
def test_gene_constraint_exposes_missense_columns():
    # NIPBL is strongly missense-constrained (mis_z ~5.6) → PP2 gene, not a BP1 gene.
    nipbl = extra.gene_constraint("NIPBL")
    assert nipbl["missense_constrained"] is True
    assert nipbl["mis_z"] is not None and nipbl["mis_z"] >= extra.MIS_Z_CONSTRAINED
    assert nipbl["missense_tolerant"] is False


def test_gene_constraint_null_for_unknown_gene():
    row = extra.gene_constraint("NOT_A_REAL_GENE_XYZ")
    assert row["missense_constrained"] is None
    assert row["missense_tolerant"] is None
    assert row["mis_z"] is None
