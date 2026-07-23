"""PS1 / PM5 — residue-level ClinVar cross-match (now engine-decided).

* **PS1** (Strong): same amino-acid change as a *different* established ClinVar pathogenic
  variant (the query's own record is PP5, not PS1).
* **PM5** (Moderate): a *different* pathogenic missense at the same residue, applied only
  when the query's exact change is not itself established.

The two are mutually exclusive by construction. Both read the residue index built by
``scripts/fetch_clinvar_residue.py`` (or the committed frozen slice).
"""
import gzip

from vcf2report.acmg.criteria import all_criteria
from vcf2report.annotate import clinvar_residue
from vcf2report.models import Annotation, Variant

_ps1 = all_criteria()["PS1"]
_pm5 = all_criteria()["PM5"]


def _v(hgvs_p="p.Arg123Cys", gene="TESTG"):
    return Variant(chrom="1", pos=100, ref="C", alt="T", gene=gene,
                   consequence="missense_variant", hgvs_p=hgvs_p)


# --- criterion logic (Annotation-driven) ------------------------------------
def _ann(**kw):
    kw.setdefault("clinvar_residue_available", True)
    return Annotation(**kw)


def test_ps1_met_on_same_aa_different_variant():
    m = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-101-G-A",
         "accession": "VCV000012345"}
    cr = _ps1(_v(), _ann(clinvar_ps1=m))
    assert cr.applies and cr.met
    assert cr.applied_strength == "strong"
    assert cr.adjudicated_by == "engine"
    assert cr.citation == ["VCV000012345"]


def test_ps1_not_met_without_match():
    cr = _ps1(_v(), _ann(clinvar_ps1=None))
    assert cr.applies and not cr.met
    assert "no distinct ClinVar pathogenic" in cr.reasoning


def test_ps1_unavailable_index_is_honest():
    cr = _ps1(_v(), Annotation(clinvar_residue_available=False))
    assert not cr.met
    assert "unavailable" in cr.reasoning
    assert cr.confidence == "low"


def test_pm5_met_on_different_aa_same_residue():
    m = {"alt_aa": "His", "ref_aa": "Arg", "stars": 1, "genomic_key": "1-101-G-A",
         "accession": "VCV000067890"}
    cr = _pm5(_v(), _ann(clinvar_pm5=m, clinvar_ps1=None))
    assert cr.applies and cr.met
    assert cr.applied_strength == "moderate"
    assert cr.citation == ["VCV000067890"]


def test_pm5_suppressed_when_ps1_fires():
    # If the exact change is itself established (PS1), PM5 must not also fire.
    ps1 = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-101-G-A"}
    pm5 = {"alt_aa": "His", "ref_aa": "Arg", "stars": 1, "genomic_key": "1-102-G-A"}
    cr = _pm5(_v(), _ann(clinvar_ps1=ps1, clinvar_pm5=pm5))
    assert not cr.met
    assert "PS1" in cr.reasoning


def test_ps1_pm5_not_applicable_for_non_missense():
    v = _v(hgvs_p=None)
    assert not _ps1(v, _ann()).met
    assert not _pm5(v, _ann()).met


# --- loader + lookup (real table format) ------------------------------------
def _write_index(path, rows):
    with gzip.open(path, "wt") as w:
        w.write("# ClinVar residue index for PS1/PM5\n")
        w.write("# Columns: gene\taa_pos\tref_aa\talt_aa\tstars\tgenomic_key\taccession\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")


def _reset_index(monkeypatch, frozen=None):
    from vcf2report import config
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_FROZEN", frozen)
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_LOCAL", None)
    clinvar_residue._index = None


def test_lookup_ps1_and_pm5_from_table(tmp_path, monkeypatch):
    fp = tmp_path / "residue.tsv.gz"
    _write_index(fp, [
        # BRCA1 Arg1699: two distinct pathogenic AA changes at the same residue.
        ("BRCA1", 1699, "Arg", "Trp", 3, "17-43057051-G-A", "VCV000000001"),
        ("BRCA1", 1699, "Arg", "Gln", 2, "17-43057052-C-T", "VCV000000002"),
    ])
    _reset_index(monkeypatch, frozen=fp)

    # Query = same AA change (Arg1699Trp) via a DIFFERENT nucleotide -> PS1.
    r = clinvar_residue.lookup("BRCA1", "p.Arg1699Trp", "17-99999999-A-T")
    assert r["available"] and r["ps1"] and r["ps1"]["accession"] == "VCV000000001"
    assert r["pm5"] is None  # exact change known -> not PM5

    # Query = novel AA change (Arg1699Leu) at the same residue -> PM5 (strongest = Trp, 3★).
    r = clinvar_residue.lookup("BRCA1", "p.Arg1699Leu", "17-99999999-A-T")
    assert r["ps1"] is None and r["pm5"] and r["pm5"]["stars"] == 3


def test_lookup_query_own_record_is_not_ps1(tmp_path, monkeypatch):
    fp = tmp_path / "residue.tsv.gz"
    _write_index(fp, [("SCN1A", 123, "Arg", "Cys", 2, "2-165991000-G-A", "VCV000000009")])
    _reset_index(monkeypatch, frozen=fp)
    # Same genomic key as the stored record -> the variant's OWN assertion (PP5), never PS1.
    r = clinvar_residue.lookup("SCN1A", "p.Arg123Cys", "2-165991000-G-A")
    assert r["ps1"] is None and r["pm5"] is None


def test_lookup_unavailable_when_no_table(monkeypatch):
    _reset_index(monkeypatch, frozen=None)
    r = clinvar_residue.lookup("BRCA1", "p.Arg1699Trp", "17-1-A-T")
    assert r["available"] is False and r["ps1"] is None and r["pm5"] is None


def test_parse_hgvs_p_rejects_non_missense():
    assert clinvar_residue.parse_hgvs_p("p.Arg123Cys") == ("Arg", 123, "Cys")
    assert clinvar_residue.parse_hgvs_p("p.Arg123Ter") is None   # nonsense
    assert clinvar_residue.parse_hgvs_p("p.Leu479fs") is None     # frameshift
    assert clinvar_residue.parse_hgvs_p("p.Ser330Ser") is None    # synonymous
    assert clinvar_residue.parse_hgvs_p(None) is None
