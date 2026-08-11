"""Gene->disease store: per-DISEASE inheritance (build + load), mocked HPO inputs."""
import subprocess
import sys
from pathlib import Path

from vcf2report import config
from vcf2report.annotate import gene_disease

ROOT = Path(__file__).resolve().parent.parent

# The real-world case from the improvement notes: JUP causes two diseases with DIFFERENT
# inheritance. A per-gene label would collapse them to "AD/AR" and hide what matters.
G2P = """\
ncbi_gene_id\tgene_symbol\thpo_id\thpo_name\tfrequency\tdisease_id
3728\tJUP\tHP:0000007\tAutosomal recessive inheritance\t-\tOMIM:601214
3728\tJUP\tHP:0002617\tDilatation\t-\tOMIM:601214
3728\tJUP\tHP:0000006\tAutosomal dominant inheritance\t-\tOMIM:611528
3728\tJUP\tHP:0001645\tSudden cardiac death\t-\tOMIM:611528
"""
HPOA = """\
#description: test
database_id\tdisease_name\tqualifier\thpo_id\treference\tevidence\tonset\tfrequency\tsex\tmodifier\taspect\tbiocuration
OMIM:601214\tNaxos disease\t\tHP:0000007\t\t\t\t\t\t\t\t
OMIM:611528\tArrhythmogenic right ventricular dysplasia 12\t\tHP:0000006\t\t\t\t\t\t\t\t
"""


def test_build_and_load_inheritance_is_per_disease(tmp_path, monkeypatch):
    g2p = tmp_path / "genes_to_phenotype.txt"
    hpoa = tmp_path / "phenotype.hpoa"
    store = tmp_path / "gene_disease.json.gz"
    g2p.write_text(G2P)
    hpoa.write_text(HPOA)

    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_gene_disease_store.py"),
         str(g2p), str(hpoa), "--out", str(store)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert store.exists()

    # load through the runtime path
    monkeypatch.setattr(config, "HPO_GENE_DISEASE_LOCAL", store)
    gene_disease._store = None  # reset the module cache
    diseases = {d["disease_name"]: d["inheritance"] for d in gene_disease.lookup("JUP")}

    assert diseases["Naxos disease"] == ["Autosomal recessive"]
    assert diseases["Arrhythmogenic right ventricular dysplasia 12"] == ["Autosomal dominant"]
    # case-insensitive gene lookup; unknown gene -> empty (never a false "no disease")
    assert gene_disease.lookup("jup")
    assert gene_disease.lookup("NOTAGENE") == []


def test_lookup_empty_when_store_absent(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "HPO_GENE_DISEASE_LOCAL", tmp_path / "missing.json.gz")
    gene_disease._store = None
    assert gene_disease.lookup("JUP") == []
