"""CLI HPO-term resolution: a file OR an inline HP: list — never silently dropped
(an inline list treated as a bogus path would disable PP4 without warning)."""
from vcf2report.cli import read_hpo_file


def test_reads_hpo_from_file(tmp_path):
    f = tmp_path / "hpo.txt"
    f.write_text("HP:0001250  seizure\n# comment\n\nHP:0001263\n")
    assert read_hpo_file(str(f)) == ["HP:0001250", "HP:0001263"]


def test_parses_inline_comma_list():
    assert read_hpo_file("HP:0001250,HP:0001263") == ["HP:0001250", "HP:0001263"]


def test_parses_inline_space_and_lowercase():
    assert read_hpo_file("hp:0001250 HP:0001263") == ["HP:0001250", "HP:0001263"]


def test_missing_plain_path_is_empty():
    assert read_hpo_file("/no/such/file.txt") == []
