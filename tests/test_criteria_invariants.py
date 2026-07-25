"""Cross-criterion invariants — the rules that must hold no matter what the inputs are.

ACMG's combining rules assume each line of evidence is counted ONCE. Two double-counts have already
slipped in through otherwise-correct individual criteria:

* a variant's own ClinVar assertion earning **PP5 and PS1** (same ClinVar record, twice);
* the same missense-intolerance signal earning **PP2 and PM1** (gene-wide and region, twice).

Both were invisible in the per-criterion tests, because each criterion was individually right. These
sweeps assert the *relationships* instead, across value ranges rather than single fixtures, so a
future change to any one criterion cannot silently re-open a gap.
"""
import pytest

from vcf2report.acmg.engine import evaluate_criteria
from vcf2report.models import Annotation, Variant


def _fired(v, a) -> set[str]:
    return {c.code for c in evaluate_criteria(v, a) if c.met}


def _missense(gene="TESTG", hgvs_p="p.Arg123Cys"):
    return Variant(chrom="1", pos=100, ref="C", alt="T", gene=gene,
                   consequence="missense_variant", hgvs_p=hgvs_p)


_HOTSPOT = {"n_residues": 5, "n_changes": 12, "window": 7, "enrichment": 4.0,
            "gene_baseline": 0.05, "available": True}
_PM5 = {"alt_aa": "His", "ref_aa": "Arg", "stars": 2, "n_other": 1, "accession": "VCV2"}
_PS1 = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-999-A-T",
        "accession": "VCV1"}


# --- in-silico: PP3 and BP4 point opposite ways -----------------------------
@pytest.mark.parametrize("am", [i / 100 for i in range(0, 101, 5)])
def test_pp3_and_bp4_never_both_fire_on_alphamissense(am):
    assert not {"PP3", "BP4"} <= _fired(_missense(), Annotation(am_pathogenicity=am))


@pytest.mark.parametrize("revel", [0.0, 0.1, 0.15, 0.5, 0.7, 0.9, 1.0])
@pytest.mark.parametrize("cadd", [0, 5, 10, 15, 20, 30])
def test_pp3_and_bp4_never_both_fire_on_revel_cadd(revel, cadd):
    fired = _fired(_missense(), Annotation(revel=revel, cadd_phred=cadd))
    assert not {"PP3", "BP4"} <= fired


# --- frequency: rarity and commonness are exclusive -------------------------
@pytest.mark.parametrize("af", [0.0, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.03, 0.05, 0.1, 0.5])
def test_at_most_one_frequency_criterion_fires(af):
    fired = _fired(_missense(), Annotation(gnomad_af=af, gnomad_faf95=af, abraom_af=0.0))
    assert len(fired & {"PM2", "BS1", "BA1"}) <= 1, f"af={af} fired {fired}"


# --- null vs protein-length: PVS1 and PM4 partition the consequences --------
@pytest.mark.parametrize("consequence", [
    "stop_gained", "frameshift_variant", "splice_donor_variant", "splice_acceptor_variant",
    "start_lost", "stop_lost", "inframe_deletion", "inframe_insertion",
    "disruptive_inframe_deletion", "conservative_inframe_insertion",
])
def test_pvs1_and_pm4_never_both_fire(consequence):
    v = Variant(chrom="1", pos=100, ref="C", alt="T", gene="TESTG",
                consequence=consequence, exon="3/10")
    fired = _fired(v, Annotation(gene_lof_intolerant=True))
    assert not {"PVS1", "PM4"} <= fired, f"{consequence} fired {fired}"


# --- the ClinVar residue hierarchy: same change / same residue / neighbourhood ---
@pytest.mark.parametrize("ann_kw, expected", [
    ({"clinvar_ps1": _PS1}, "PS1"),                      # same residue, same AA change
    ({"clinvar_pm5": _PM5}, "PM5"),                      # same residue, different AA change
    ({"clinvar_hotspot": _HOTSPOT}, "PM1"),              # neighbouring residues only
])
def test_residue_evidence_fires_exactly_one_criterion(ann_kw, expected):
    a = Annotation(clinvar_residue_available=True, **ann_kw)
    fired = _fired(_missense(), a) & {"PS1", "PM5", "PM1"}
    assert fired == {expected}, f"{ann_kw} fired {fired}"


def test_residue_hierarchy_holds_when_every_signal_is_present():
    # All three signals at once: PS1 (most specific) wins, PM5 and PM1 both stand down, so the same
    # residue evidence is never counted two or three times.
    a = Annotation(clinvar_residue_available=True, clinvar_ps1=_PS1, clinvar_pm5=_PM5,
                   clinvar_hotspot=_HOTSPOT)
    assert _fired(_missense(), a) & {"PS1", "PM5", "PM1"} == {"PS1"}


# --- gene-level vs region-level missense intolerance ------------------------
def test_pp2_and_pm1_never_both_fire():
    a = Annotation(clinvar_residue_available=True, clinvar_hotspot=_HOTSPOT,
                   gene_mis_z=5.5, gene_missense_constrained=True)
    fired = _fired(_missense(), a)
    assert "PM1" in fired and "PP2" not in fired


def test_pp2_and_bp1_never_both_fire():
    # A gene cannot be both depleted of and tolerant to missense; assert it even if flags disagree.
    a = Annotation(gene_missense_constrained=True, gene_missense_tolerant=True,
                   gene_lof_intolerant=True, gene_mis_z=5.5)
    assert not {"PP2", "BP1"} <= _fired(_missense(), a)


# --- a variant's own ClinVar record is counted once --------------------------
@pytest.mark.parametrize("significance", ["Pathogenic", "Likely pathogenic"])
def test_own_clinvar_assertion_never_earns_pp5_and_a_residue_criterion(significance):
    a = Annotation(clinvar_residue_available=True, clinvar_ps1=_PS1, clinvar_pm5=_PM5,
                   clinvar_significance=significance,
                   clinvar_review_status="criteria provided, single submitter")
    fired = _fired(_missense(), a)
    assert "PP5" in fired
    assert not (fired & {"PS1", "PM5"}), f"{significance} double-counted: {fired}"


def test_pp5_and_bp6_never_both_fire():
    for sig in ("Pathogenic", "Benign", "Likely benign", "Likely pathogenic"):
        a = Annotation(clinvar_significance=sig,
                       clinvar_review_status="criteria provided, single submitter")
        assert not {"PP5", "BP6"} <= _fired(_missense(), a), sig


# --- nothing fires on absent evidence ---------------------------------------
def test_empty_annotation_fires_no_pathogenic_criterion():
    # An annotation with nothing looked up must not manufacture evidence — the honesty invariant
    # (PM2 in particular must not read "absent" out of a failed lookup).
    fired = _fired(_missense(), Annotation())
    assert not fired, f"criteria fired on an empty annotation: {fired}"
