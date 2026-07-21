"""The realistic (tell-free) spike builder: a planted record must be indistinguishable from a real
background call — no SPIKED/GENE/CSQ/CLN* markers, the full DRAGEN FORMAT preserved, the genotype
matching the requested zygosity — while still landing at the spike's own coordinate/alleles."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spike_variant import spiked_line_realistic  # noqa: E402

_FMT = "GT:AD:AF:DP:F1R2:F2R1:GQ:PL:GP:PRI:SB:MB"
# Two real-shaped DRAGEN background calls to borrow from: one het, one hom.
_HET = ["chr1", "111", ".", "G", "A", "54.1", "PASS",
        "AC=1;AF=0.5;AN=2;DP=42;FS=1.2;MQ=250;QD=1.5;SOR=0.6;ANN=A|missense_variant|MODERATE|FOO|x",
        _FMT, "0/1:23,19:0.45:42:11,10:12,9:62:98,0,71:1,2,3:0,34,37:10,13,10,9:14,9,5,14"]
_HOM = ["chr2", "222", ".", "C", "T", "66.9", "PASS",
        "AC=2;AF=1;AN=2;DP=40;FS=0;MQ=250;QD=1.6;SOR=1.1;ANN=T|synonymous_variant|LOW|BAR|y",
        _FMT, "1/1:0,40:1:40:0,22:0,18:60:105,64,0:1,2,3:0,34,37:0,0,16,24:0,0,20,20"]
_TELLS = ("SPIKED", "GENE=", "CSQ=", "CLNSIG", "CLNREVSTAT", "CLNDN", "CLNVID")


def _rec():
    return {"chrom": "7", "pos": 900900, "ref": "A", "alt": "G", "gene": "TARGET"}


def test_realistic_het_is_tell_free_and_well_formed():
    out = spiked_line_realistic(_rec(), [_HET, _HOM], "chr", "het", 10)
    assert out is not None
    assert out[:5] == ["chr7", "900900", ".", "A", "G"]          # spike's own locus/alleles
    assert out[6] == "PASS"
    assert not any(t in out[7] for t in _TELLS)                   # no markers
    assert "ANN=" not in out[7] and "LOF=" not in out[7]         # borrowed annotation stripped
    assert out[8] == _FMT and "F1R2" in out[8]                    # full DRAGEN FORMAT kept
    assert out[9].split(":")[0] == "0/1"                          # het genotype


def test_realistic_hom_borrows_a_hom_template():
    out = spiked_line_realistic(_rec(), [_HET, _HOM], "chr", "hom", 10)
    assert out is not None
    assert out[9].split(":")[0] == "1/1"
    assert "AC=2" in out[7] and "AN=2" in out[7]                  # hom call stats, not het
    assert not any(t in out[7] for t in _TELLS)


def test_realistic_returns_none_without_a_matching_template():
    # No hom template present → cannot borrow a realistic hom call.
    assert spiked_line_realistic(_rec(), [_HET], "chr", "hom", 10) is None
    # A minimal-FORMAT record (no F1R2) is not a usable template.
    minimal = ["chr1", "5", ".", "G", "A", "800", "PASS", "GENE=X;SPIKED=1", "GT:DP:GQ:AD", "0/1:44:99:22,22"]
    assert spiked_line_realistic(_rec(), [minimal], "chr", "het", 10) is None


def test_plain_style_has_no_chr_prefix():
    out = spiked_line_realistic(_rec(), [_HET, _HOM], "plain", "het", 10)
    assert out[0] == "7"
