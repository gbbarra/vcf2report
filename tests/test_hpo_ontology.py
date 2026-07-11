"""Ontology-aware HPO matching (Lin/IC) with exact-overlap fallback."""
import pytest

from vcf2report import config
from vcf2report.annotate import hpo

# Tiny ontology: Seizure(1.0) -> {Focal(2.0), Status(2.5)}; Cancer(1.2). Root IC 0.
GRAPH = (
    "# id\tname\tic\tparents\n"
    "HP:0\tRoot\t0\t\n"
    "HP:1\tSeizure\t1.0\tHP:0\n"
    "HP:11\tFocal seizure\t2.0\tHP:1\n"
    "HP:12\tStatus epilepticus\t2.5\tHP:1\n"
    "HP:2\tCancer\t1.2\tHP:0\n"
)


@pytest.fixture
def graph(tmp_path, monkeypatch):
    fp = tmp_path / "g.tsv"
    fp.write_text(GRAPH)
    monkeypatch.setattr(config, "HPO_GRAPH_LOCAL", fp)
    monkeypatch.setattr(hpo, "_graph", None)
    yield
    monkeypatch.setattr(hpo, "_graph", None)


def test_lin_similarity(graph):
    parents, ic, names = hpo._load_graph()
    cache = {}
    L = lambda a, b: round(hpo._lin(a, b, parents, ic, cache), 3)
    assert L("HP:11", "HP:11") == 1.0                     # identity
    assert L("HP:1", "HP:11") == 0.667                    # parent<->child
    assert L("HP:11", "HP:12") == 0.444                   # siblings share Seizure
    assert L("HP:11", "HP:2") == 0.0                      # only the root in common


def test_bma_no_dilution(graph):
    # Two patient terms both exactly in the gene spectrum -> perfect.
    assert hpo._semantic_match({"HP:11", "HP:12"}, ["HP:11", "HP:12"])["score"] == 1.0
    # A third patient term that is a (related) parent barely moves the score — the
    # exact-overlap failure mode (more HPO -> lower score) is gone.
    r = hpo._semantic_match({"HP:11", "HP:12"}, ["HP:11", "HP:12", "HP:1"])
    assert r["score"] >= 0.8
    # An unrelated gene (cancer) vs a seizure patient stays ~0.
    assert hpo._semantic_match({"HP:2"}, ["HP:11", "HP:12"])["score"] == 0.0


def test_related_credit_beats_exact_overlap(graph):
    # Gene annotated only with Focal seizure; patient has Status epilepticus (sibling).
    # Exact overlap would score 0; the ontology credits the related term.
    r = hpo._semantic_match({"HP:11"}, ["HP:12"])
    assert 0 < r["score"] < 1
    assert "HP:12→HP:11" in r["matched_terms"][0]


def test_fallback_exact_overlap_when_no_graph(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "HPO_GRAPH_LOCAL", tmp_path / "absent.tsv")
    genes = tmp_path / "genes.tsv"
    genes.write_text("G\tHP:11\tFocal seizure\nG\tHP:99\tOther\n")
    monkeypatch.setattr(config, "HPO_GENES_LOCAL", genes)
    monkeypatch.setattr(hpo, "_graph", None)
    monkeypatch.setattr(hpo, "_gene_terms", None)
    monkeypatch.setattr(hpo, "_term_names", {})
    r = hpo.match("G", ["HP:11", "HP:88"])   # 1 of 2 patient terms overlaps exactly
    assert r["_source"].startswith("HPO genes_to_phenotype")
    assert r["score"] == 0.5


def test_real_graph_discriminates(monkeypatch):
    if not config.HPO_GRAPH_LOCAL.exists():
        pytest.skip("HPO graph not built")
    monkeypatch.setattr(hpo, "_graph", None)
    monkeypatch.setattr(hpo, "_gene_terms", None)
    epilepsy = ["HP:0001250", "HP:0002133", "HP:0011097"]
    assert hpo.match("SCN1A", epilepsy)["score"] >= 0.9          # clearly related
    assert hpo.match("BRCA1", epilepsy)["score"] < config.HPO_RELATED_MIN  # incidental
