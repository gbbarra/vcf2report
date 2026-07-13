"""The full local ClinVar store (tabix) is the authoritative offline source: a real assertion
is returned without any network call; a miss falls through to 'no record' (offline)."""
import pytest

pysam = pytest.importorskip("pysam")

from vcf2report import config
from vcf2report.annotate import clinvar
from vcf2report.models import Variant


def _store(tmp_path):
    raw = tmp_path / "cv.tsv"
    raw.write_text(
        "1\t100\tA\tT\tPathogenic\tcriteria provided, multiple submitters, no conflicts\t12345\tSome disease\n"
        "1\t200\tC\tG\tBenign\tcriteria provided, single submitter\t67890\tOther condition\n")
    gz = str(raw) + ".gz"
    pysam.tabix_compress(str(raw), gz, force=True)
    pysam.tabix_index(gz, seq_col=0, start_col=1, end_col=1, force=True)
    return gz


def test_local_clinvar_is_authoritative_offline(tmp_path, monkeypatch):
    monkeypatch.setenv("OFFLINE", "1")
    monkeypatch.setattr(config, "CLINVAR_TABIX", _store(tmp_path))
    clinvar._tabix = None
    try:
        r = clinvar.lookup(Variant(chrom="chr1", pos=100, ref="A", alt="T", gene="X", consequence="missense_variant"))
        assert r["significance"] == "Pathogenic"
        assert "multiple submitters" in r["review_status"]
        assert r["_source"] == "ClinVar (local)"
        # a locus not in the store -> no record (offline, and the tiny slice doesn't cover it)
        miss = clinvar.lookup(Variant(chrom="chr1", pos=999, ref="A", alt="T", gene="X", consequence="missense_variant"))
        assert miss["significance"] is None
    finally:
        clinvar._tabix = None
