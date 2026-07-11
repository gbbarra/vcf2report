"""Safety guard: a configured-but-unavailable gnomAD parquet must be flagged loudly.

If the operator points VCF2REPORT_GNOMAD_PARQUET at a store that is missing (e.g. an
unmounted drive), silent fall-through would make every variant look absent from gnomAD
-> spurious PM2 everywhere -> gross over-calling. The pipeline must warn instead.
"""
import pytest

from vcf2report import config, pipeline

pytest.importorskip("duckdb")

_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##reference=GRCh38\n"
    '##INFO=<ID=GENE,Number=1,Type=String,Description="g">\n'
    '##INFO=<ID=CSQ,Number=1,Type=String,Description="c">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">\n'
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="dp">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="gq">\n'
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="ad">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def _big_vcf(tmp_path, n=60):
    rows = [_HEADER]
    for i in range(n):
        pos = 100000 + i * 137
        rows.append(
            f"1\t{pos}\t.\tA\tG\t800\tPASS\tGENE=BRCA1;CSQ=missense_variant\t"
            "GT:DP:GQ:AD\t0/1:40:99:20,20\n"
        )
    p = tmp_path / "big.vcf"
    p.write_text("".join(rows))
    return str(p)


def _reset():
    from vcf2report.annotate import gnomad_parquet
    gnomad_parquet._reset_for_tests()


def test_configured_but_missing_parquet_warns(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(tmp_path / "nope.parquet"))
    _reset()
    report = pipeline.run_pipeline(_big_vcf(tmp_path))
    assert any("resolved from it" in w and "OVER-calls" in w for w in report.qc.warnings)
    _reset()


def test_unconfigured_parquet_does_not_warn(tmp_path, monkeypatch):
    monkeypatch.setattr(config, "GNOMAD_PARQUET", None)
    _reset()
    report = pipeline.run_pipeline(_big_vcf(tmp_path))
    assert not any("OVER-calls" in w for w in report.qc.warnings)
    _reset()


def test_present_store_but_duckdb_missing_blames_duckdb(tmp_path, monkeypatch):
    # Store dir EXISTS (e.g. fetched with gh+zstd) but duckdb isn't installed -> the
    # guard must point at `pip install duckdb`, not an unmounted drive.
    store = tmp_path / "gnomad_parquet"
    store.mkdir()
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(store))
    from vcf2report.annotate import gnomad_parquet
    monkeypatch.setattr(gnomad_parquet, "_get_duckdb", lambda: None)
    _reset()
    report = pipeline.run_pipeline(_big_vcf(tmp_path))
    assert any("pip install duckdb" in w for w in report.qc.warnings)
    assert not any("unmounted" in w for w in report.qc.warnings)
    _reset()


def test_small_callset_does_not_warn(tmp_path, monkeypatch):
    # A tiny demo VCF (< 50 kept) must not trip the guard even if the store is absent.
    monkeypatch.setattr(config, "GNOMAD_PARQUET", str(tmp_path / "nope.parquet"))
    _reset()
    report = pipeline.run_pipeline(_big_vcf(tmp_path, n=10))
    assert not any("OVER-calls" in w for w in report.qc.warnings)
    _reset()
