"""Structural (network-free) tests for scripts/build_gnomad_parquet.py.

The per-population AF columns are keyed to exact gnomAD INFO field names; a wrong
name would not crash — bcftools emits '.' and the column silently becomes all-NULL.
These tests pin the field mapping so that regression is caught without a network build.
"""
import importlib.util
import json
import pathlib

import pytest

_SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "build_gnomad_parquet.py"
_spec = importlib.util.spec_from_file_location("build_gnomad_parquet", _SCRIPT)
bgp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bgp)


def test_resume_skips_only_on_scope_match(tmp_path):
    # A partition is reused ONLY if its _scope.json matches; else a later broader build
    # (dropping --bed / switching preset) must re-extract, never silently reuse a slice.
    duckdb = pytest.importorskip("duckdb")
    out = tmp_path / "store"
    part = out / "chrom=chr21"
    part.mkdir(parents=True)
    duckdb.connect().execute(
        f"COPY (SELECT 1 AS chrom, 100 AS pos) TO '{part/'data.parquet'}' (FORMAT PARQUET)")
    (part / "_scope.json").write_text(json.dumps({"region": None, "bed": None, "preset": "exomes"}))
    url, fields = bgp._PRESETS["exomes"]

    # same scope -> resume-skip, returns the existing row count (never touches the source)
    assert bgp.build_chrom("21", url, fields, "/nonexistent", out, None, None, None, "exomes") == 1
    # different scope (preset) -> must NOT skip; falls through to a missing source -> 0
    assert bgp.build_chrom("21", url, fields, "/nonexistent", out, None, None, None, "joint") == 0


def _fields(preset):
    return dict(bgp._PRESETS[preset][1])


def test_joint_has_all_ten_populations():
    f = _fields("joint")
    for pop in ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]:
        assert f[f"af_{pop}"] == f"AF_joint_{pop}"


def test_exomes_populations_exclude_ami():
    f = _fields("exomes")
    assert "af_ami" not in f            # the exomes-only release has no Amish group
    for pop in ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "remaining", "sas"]:
        assert f[f"af_{pop}"] == f"AF_{pop}"


def test_core_columns_still_present():
    for preset in ("joint", "exomes"):
        f = _fields(preset)
        assert {"af", "af_grpmax", "ac", "an", "nhomalt", "faf95", "grpmax_pop"} <= set(f)


def test_every_pop_af_is_numeric_double():
    for preset in ("joint", "exomes"):
        for out, _info in bgp._PRESETS[preset][1]:
            if out.startswith("af_"):
                assert bgp._NUMERIC.get(out) == "DOUBLE"


def test_bcftools_fmt_and_copy_sql_reference_the_pop_fields():
    fields = bgp._PRESETS["joint"][1]
    fmt = bgp._bcftools_fmt(fields)
    assert "%INFO/AF_joint_nfe" in fmt and "%INFO/AF_joint_afr" in fmt
    cols = bgp._LEAD + [out for out, _ in fields]
    sql = bgp._copy_sql(pathlib.Path("/tmp/x.tsv"), pathlib.Path("/tmp/o.parquet"), cols)
    # per-pop AFs must be cast to DOUBLE (not left as VARCHAR)
    assert "AS af_nfe" in sql and "TRY_CAST(NULLIF(af_nfe" in sql
