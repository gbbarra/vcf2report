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


def test_best_single_is_undiluted(graph):
    # A gene that explains ONE of several patient terms perfectly: the average is
    # diluted (drives PP4) but 'best' stays 1.0 (drives primary routing).
    r = hpo._semantic_match({"HP:11"}, ["HP:11", "HP:2"])
    assert r["best"] == 1.0
    assert r["score"] == 0.5


def test_fallback_best_on_any_overlap(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "HPO_GRAPH_LOCAL", tmp_path / "absent.tsv")
    genes = tmp_path / "genes.tsv"
    genes.write_text("G\tHP:11\tFocal seizure\n")
    monkeypatch.setattr(config, "HPO_GENES_LOCAL", genes)
    monkeypatch.setattr(hpo, "_graph", None)
    monkeypatch.setattr(hpo, "_gene_terms", None)
    monkeypatch.setattr(hpo, "_term_names", {})
    # A single exact overlap keeps best=1.0 so routing still puts it in primary.
    r = hpo.match("G", ["HP:11", "HP:99", "HP:98"])
    assert r["score"] < 0.5 and r["best"] == 1.0


def test_ancestors_cycle_safe():
    # A corrupted cyclic graph must degrade, not overflow the stack.
    parents = {"A": ["B"], "B": ["A"], "C": ["C"]}
    assert hpo._ancestors("A", parents, {}) == {"A", "B"}
    assert hpo._ancestors("C", parents, {}) == {"C"}


def test_defect1_unannotated_ic_is_high_not_zero():
    import importlib.util
    from pathlib import Path
    fp = Path(__file__).resolve().parent.parent / "scripts" / "build_hpo_graph.py"
    spec = importlib.util.spec_from_file_location("build_hpo_graph", fp)
    bhg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bhg)
    parents = {"HP:0": [], "HP:1": ["HP:0"], "HP:11": ["HP:1"], "HP:12": ["HP:1"], "HP:2": ["HP:0"]}
    ic = bhg.compute_ic(parents, [("G1", "HP:11"), ("G2", "HP:2")])
    assert ic["HP:12"] > 0                  # unannotated term is specific, not IC=0
    assert ic["HP:12"] >= ic["HP:1"]        # IC monotonic non-decreasing downward
    assert ic["HP:0"] == 0.0                # root stays uninformative


def test_real_graph_discriminates(monkeypatch):
    if not config.HPO_GRAPH_LOCAL.exists():
        pytest.skip("HPO graph not built")
    monkeypatch.setattr(hpo, "_graph", None)
    monkeypatch.setattr(hpo, "_gene_terms", None)
    epilepsy = ["HP:0001250", "HP:0002133", "HP:0011097"]
    assert hpo.match("SCN1A", epilepsy)["score"] >= 0.9          # clearly related
    assert hpo.match("BRCA1", epilepsy)["score"] < config.HPO_RELATED_MIN  # incidental
