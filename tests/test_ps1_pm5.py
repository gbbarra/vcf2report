"""PS1 / PM5 — residue-level ClinVar cross-match (now engine-decided).

* **PS1** (Strong): same amino-acid change as a *different* established ClinVar pathogenic
  variant (the query's own record is PP5, not PS1).
* **PM5** (Moderate): a *different* pathogenic missense at the same residue, applied only
  when the query's exact change is not itself established.

The two are mutually exclusive by construction. Both read the residue index built by
``scripts/fetch_clinvar_residue.py`` (or the committed frozen slice).
"""
import pytest

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
    # One other pathogenic change, well-reviewed (>=2*) -> PM5 at its default Moderate.
    # (The graded Strong/Supporting variants are covered in tests/test_pm1_hotspot.py.)
    m = {"alt_aa": "His", "ref_aa": "Arg", "stars": 2, "n_other": 1,
         "genomic_key": "1-101-G-A", "accession": "VCV000067890"}
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


def test_ps1_suppressed_when_own_clinvar_is_pathogenic():
    # A variant ClinVar already calls Pathogenic is covered by PP5; PS1 (which would also fire
    # from a same-AA different-locus entry) is withheld so the ClinVar evidence isn't double-counted.
    m = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-101-G-A",
         "accession": "VCV000012345"}
    cr = _ps1(_v(), _ann(clinvar_ps1=m, clinvar_significance="Pathogenic",
                         clinvar_review_status="criteria provided, single submitter"))
    assert not cr.met
    assert "PP5" in cr.reasoning and "double-count" in cr.reasoning


def test_pm5_suppressed_when_own_clinvar_is_pathogenic():
    m = {"alt_aa": "His", "ref_aa": "Arg", "stars": 1, "genomic_key": "1-101-G-A",
         "accession": "VCV000067890"}
    cr = _pm5(_v(), _ann(clinvar_pm5=m, clinvar_ps1=None,
                         clinvar_significance="Likely pathogenic",
                         clinvar_review_status="criteria provided, single submitter"))
    assert not cr.met
    assert "PP5" in cr.reasoning


def test_ps1_not_suppressed_by_an_unreviewed_clinvar_record():
    # The guard's justification is "captured by PP5" — so it must test the condition PP5 actually
    # tests, stars included. A 0-star "Pathogenic" does NOT fire PP5, so withholding PS1 on it
    # would strip a legitimate residue cross-match while citing a PP5 that never fired.
    m = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-101-G-A",
         "accession": "VCV000012345"}
    a = _ann(clinvar_ps1=m, clinvar_significance="Pathogenic",
             clinvar_review_status="no assertion criteria provided")
    assert not all_criteria()["PP5"](_v(), a).met      # 0-star: PP5 does not fire...
    assert _ps1(_v(), a).met                           # ...so PS1 must not stand down


def test_ps1_still_fires_when_own_clinvar_is_vus():
    # A variant ClinVar calls VUS (not P/LP) still earns PS1 from a same-AA pathogenic elsewhere —
    # PP5 does not fire for a VUS, so there is nothing to double-count.
    m = {"alt_aa": "Cys", "ref_aa": "Arg", "stars": 2, "genomic_key": "1-101-G-A",
         "accession": "VCV000012345"}
    cr = _ps1(_v(), _ann(clinvar_ps1=m, clinvar_significance="Uncertain significance"))
    assert cr.met


# --- loader + lookup (real table format) ------------------------------------
def _write_index(path, rows):
    with gzip.open(path, "wt") as w:
        w.write("# ClinVar residue index for PS1/PM5\n")
        w.write("# Columns: gene\taa_pos\tref_aa\talt_aa\tstars\tgenomic_key\taccession\n")
        for r in rows:
            w.write("\t".join(map(str, r)) + "\n")


def _reset_index(monkeypatch, frozen=None):
    """Point the loader at a test table AND make the swap undo itself.

    The module caches the parsed index (and the per-gene baselines derived from it). Clearing them
    by plain assignment leaked the test table into every later test in the session — monkeypatch
    restored the config paths but not the cache, so the loader saw "already loaded" and kept the
    tmp data. Going through monkeypatch means teardown restores the previous cache too.
    """
    from vcf2report import config
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_FROZEN", frozen)
    monkeypatch.setattr(config, "CLINVAR_RESIDUE_LOCAL", None)
    monkeypatch.setattr(clinvar_residue, "_index", None)
    monkeypatch.setattr(clinvar_residue, "_baselines", {})


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


def test_annotate_variant_populates_residue_evidence_by_default():
    """A single `annotate_variant` call must return a COMPLETE annotation.

    Regression: the residue lookup was moved out of `annotate_variant` into a pipeline-only
    deferral step for performance. That silently stripped PS1/PM5/PM1 from every other entry
    point — notably the MCP `classify_variant` tool, which calls `annotate_variant` directly and
    then reported "index unavailable — build it" for an index that was present and would have
    fired. The deferral is now opt-in (`with_clinvar_residue=False`), as it already was for
    AlphaMissense, so ad-hoc callers get complete annotations and only the pipeline defers.
    """
    from vcf2report.annotate import annotate_variant
    from vcf2report.models import Variant as V

    v = V(chrom="5", pos=112838671, ref="A", alt="T", gene="APC",
          consequence="missense_variant", hgvs_p="p.Asn1026Asp")
    ann = annotate_variant(v, [])
    assert ann.clinvar_residue_available is True, "residue index not consulted by default"
    assert ann.clinvar_hotspot is not None, "hotspot not computed by default"

    # ...and the deferred form leaves them unset, for the pipeline to fill on the candidates.
    deferred = annotate_variant(v, [], with_clinvar_residue=False)
    assert deferred.clinvar_residue_available is None
    assert deferred.clinvar_hotspot is None


def test_ps1_requires_the_reference_amino_acid_to_agree(monkeypatch):
    """The index is keyed (gene, aa_pos, alt_aa) and stores ref_aa, but nothing compared it to
    the query's. A query Ala265Glu therefore earned PS1 Strong from an index row for Gly265Glu —
    and the reasoning string printed the INDEX's ref_aa, so the trail read "same amino-acid
    change (Gly→Glu)" and was internally consistent. A reviewer could not see the mismatch."""
    from vcf2report.annotate import clinvar_residue

    idx = {"GENEX": {265: {"Glu": ("Gly", 2, "5-1-A-G", "VCV000000001")}}}
    monkeypatch.setattr(clinvar_residue, "_index", idx)
    monkeypatch.setattr(clinvar_residue, "_baselines", {"GENEX": 0.0})

    assert clinvar_residue.lookup("GENEX", "p.Gly265Glu", "5-999-A-G")["ps1"]   # ref agrees
    for wrong in ("p.Ala265Glu", "p.Trp265Glu", "p.Pro265Glu"):
        r = clinvar_residue.lookup("GENEX", wrong, "5-999-A-G")
        assert r["ps1"] is None and r["pm5"] is None, f"{wrong} matched a Gly265 row"


@pytest.mark.parametrize("hgvs_p", [
    "p.Arg97GlyfsTer26",        # VEP / ClinVar long frameshift form
    "p.Ser330AsnfsTer5",
    "p.Gly12ValfsTer3",
])
def test_vep_frameshift_notation_is_not_read_as_a_missense(hgvs_p):
    """`_P_RE.search` matched the leading `p.Arg97Gly` of a frameshift. That fired PS1/PM5 on
    frameshifts AND poisoned the index: scripts/fetch_clinvar_residue.py runs the same parser
    over ClinVar's Name column, so indels were recorded as missense rows."""
    from vcf2report.annotate.clinvar_residue import parse_hgvs_p
    assert parse_hgvs_p(hgvs_p) is None


@pytest.mark.parametrize("hgvs_p,expected", [
    ("p.Gly265Glu", ("Gly", 265, "Glu")),
    ("ENSP00000350283.3:p.Gln356Arg", ("Gln", 356, "Arg")),   # VEP writes a protein prefix
    ("p.(Gly265Glu)", ("Gly", 265, "Glu")),                   # HGVS "predicted" parentheses
])
def test_real_missense_forms_still_parse(hgvs_p, expected):
    """Anchoring the regex must not lose the forms real annotators emit."""
    from vcf2report.annotate.clinvar_residue import parse_hgvs_p
    assert parse_hgvs_p(hgvs_p) == expected


def test_the_committed_residue_index_holds_no_indels():
    """Three rows in the shipped frozen slice are indels recorded as missense (a genomic key
    whose REF and ALT differ in length cannot be a substitution). They reached PS1 as Strong
    evidence. Regenerate the slice with the anchored parser if this fails."""
    import gzip
    import pathlib

    p = pathlib.Path("data/clinvar/clinvar_residue_frozen.tsv.gz")
    if not p.exists():
        pytest.skip("frozen residue slice not present")
    bad = []
    for line in gzip.open(p, "rt"):
        if line.startswith("#"):
            continue
        for tok in line.rstrip("\n").split("\t"):
            parts = tok.split("-")
            if len(parts) == 4 and parts[1].isdigit() and len(parts[2]) != len(parts[3]):
                bad.append(tok)
    assert not bad, f"non-substitution genomic keys in the residue index: {bad}"
