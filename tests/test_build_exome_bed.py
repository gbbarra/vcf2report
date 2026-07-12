"""Structural (network-free) tests for scripts/build_exome_bed.py.

The MANE filter matches the tag as a substring of a comma-list (tag=...,MANE_Select,...);
a `tag=MANE_Select` prefix match silently yields an EMPTY BED. Pin it.
"""
import importlib.util
import pathlib

_S = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "build_exome_bed.py"
_spec = importlib.util.spec_from_file_location("build_exome_bed", _S)
beb = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(beb)

_GFF = (
    "##gff-version 3\n"
    "chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_type=protein_coding;tag=basic,MANE_Select,CCDS\n"
    "chr1\tHAVANA\texon\t150\t260\t.\t+\t.\tgene_type=protein_coding;tag=basic\n"
    "chr1\tHAVANA\texon\t500\t600\t.\t+\t.\tgene_type=lncRNA;tag=MANE_Plus_Clinical\n"
    "chr1\tHAVANA\tCDS\t100\t200\t.\t+\t.\tgene_type=protein_coding;tag=MANE_Select\n"
    "chrZ\tHAVANA\texon\t1\t50\t.\t+\t.\ttag=MANE_Select\n"
)


def _gff(tmp_path):
    p = tmp_path / "g.gff3"
    p.write_text(_GFF)
    return str(p)


def test_mane_selects_only_mane_tagged(tmp_path):
    iv = list(beb.exon_intervals(_gff(tmp_path), "mane"))
    assert ("chr1", 99, 200) in iv          # MANE_Select
    assert ("chr1", 499, 600) in iv         # MANE_Plus_Clinical
    assert ("chr1", 149, 260) not in iv     # not MANE
    assert not any(c == "chrZ" for c, _s, _e in iv)   # non-standard contig dropped
    assert all(row for row in iv)           # CDS feature ignored (only 'exon')


def test_protein_coding_excludes_lncrna(tmp_path):
    iv = list(beb.exon_intervals(_gff(tmp_path), "protein_coding"))
    assert ("chr1", 499, 600) not in iv     # lncRNA exon excluded
    assert len(iv) == 2


def test_all_keeps_every_exon(tmp_path):
    assert len(list(beb.exon_intervals(_gff(tmp_path), "all"))) == 3


def test_merge_pads_and_collapses():
    assert beb.merge([("chr1", 100, 200), ("chr1", 150, 260)], pad=0) == [("chr1", 100, 260)]
    assert beb.merge([("chr1", 100, 200)], pad=50) == [("chr1", 50, 250)]
    # a gap wider than the padding stays split
    assert len(beb.merge([("chr1", 100, 200), ("chr1", 400, 500)], pad=50)) == 2
