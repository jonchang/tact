"""Tests for the placement report written by `tact_add_taxa`."""

import csv

import dendropy
import pytest

from tact import cli_add_taxa
from tact.report import MAX_SUMMARY_TAXA, PlacementReport, fmt_species_list

# G1 is not monophyletic in the backbone because G2 d falls inside it, and G4 is not sampled at all
TAXONOMY = "(((G1_a,G1_b,G1_c)G1,(G2_d,G2_e)G2)F1,((G3_f,G3_g)G3,(G4_h,G4_i)G4)F2)O;"
BACKBONE = "(((G1_a:3,G2_d:3):2,G1_b:5):5,(G3_f:4,G3_g:4):6);"


def read_rows(path):
    with open(path, encoding="utf-8") as fh:
        return {row["species"]: row for row in csv.DictReader(fh)}


@pytest.fixture
def placements(run_cli, tmp_path):
    taxonomy = tmp_path / "taxonomy.tre"
    taxonomy.write_text(TAXONOMY)
    backbone = tmp_path / "backbone.tre"
    backbone.write_text(BACKBONE)
    output = tmp_path / "out"
    result = run_cli(cli_add_taxa.main, ["--taxonomy", taxonomy, "--backbone", backbone, "--output", output, "-v"])
    return result, output


def test_rows(placements):
    _, output = placements
    rows = read_rows(f"{output}.placements.csv")
    assert rows.keys() == {"G1 c", "G2 e", "G4 h", "G4 i"}
    assert rows["G1 c"] == {
        "species": "G1 c",
        "expected_taxon": "G1",
        "placed_taxon": "F1",
        "ranks_up": "1",
        "new_clade": "",
        "skipped_taxa": "G1",
        "n_intruders": "1",
        "intruders": "G2 d",
    }
    assert rows["G2 e"]["expected_taxon"] == rows["G2 e"]["placed_taxon"] == "G2"
    assert rows["G2 e"]["skipped_taxa"] == ""
    # an unsampled genus is expected in the lowest taxon with sampled species
    for species in ("G4 h", "G4 i"):
        assert rows[species]["expected_taxon"] == rows[species]["placed_taxon"] == "F2"
        assert rows[species]["new_clade"] == "G4"
        assert rows[species]["ranks_up"] == "0"


def test_summary(placements):
    result, output = placements
    assert result.stdout.splitlines() == [
        "Added 4 species; 1 (25%) were placed in a higher taxon because the lowest taxon with species in the "
        "backbone is not monophyletic there.",
        "Most affected:",
        "  G1: 1 species placed in F1 (1 rank up; broken up by G2 d)",
        f"Placement of each added species: {output}.placements.csv",
    ]
    # the summary goes to the log file but is not repeated on stderr
    with open(f"{output}.log.txt", encoding="utf-8") as fh:
        assert "Added 4 species" in fh.read()
    assert "Added 4 species" not in result.stderr
    # the log now names the tips breaking up a taxon without -vv
    assert "G1: is not monophyletic (from G2 d)" in result.stderr


def make_report(tmp_path, taxonomy, backbone_tips):
    tree = dendropy.Tree.get(data=taxonomy, schema="newick")
    return PlacementReport(str(tmp_path / "placements.csv"), tree, backbone_tips)


def test_unlabelled_ranks_are_not_counted(tmp_path):
    with make_report(tmp_path, "((((A_a,A_b)A,(B_c)B))F)O;", {"A a", "B c"}) as report:
        report.taxon_not_monophyletic("A", ["A a", "B c"], ["A a", "A b"])
        report.placed(["A b"], "F")
    row = read_rows(tmp_path / "placements.csv")["A b"]
    assert row["ranks_up"] == "1"
    assert row["skipped_taxa"] == "A"


def test_placed_outside_expected_taxon(tmp_path):
    with make_report(tmp_path, "((A_a,A_b)A,(B_c)B)F;", {"A a", "B c"}) as report:
        with pytest.raises(ValueError, match="A b was placed in B, which does not contain its expected taxon A"):
            report.placed(["A b"], "B")


def test_intruders_are_backbone_tips(tmp_path):
    with make_report(tmp_path, "((A_a,A_b,A_c)A,(B_d,B_e)B)F;", {"A a", "A b", "B d"}) as report:
        # B e was added earlier and is not reported as an intruder
        intruders = report.taxon_not_monophyletic("A", ["A a", "A b", "B d", "B e"], ["A a", "A b", "A c"])
        assert intruders == ["B d"]


def test_summary_without_fallbacks(tmp_path):
    with make_report(tmp_path, "((A_a,A_b)A,(B_c)B)F;", {"A a", "B c"}) as report:
        report.placed(["A b"], "A")
    first_line = report.summary().splitlines()[0]
    assert first_line == "Added 1 species, all to the lowest taxon with species in the backbone."


def test_summary_without_additions(tmp_path):
    with make_report(tmp_path, "((A_a)A,(B_c)B)F;", {"A a", "B c"}) as report:
        pass
    assert report.summary().splitlines()[0] == "No species were added."


def test_summary_lists_most_affected_taxa(tmp_path):
    n_genera = MAX_SUMMARY_TAXA + 2
    genera = ",".join(f"(G{i}_a,G{i}_b)G{i}" for i in range(n_genera))
    backbone_tips = {f"G{i} a" for i in range(n_genera)}
    with make_report(tmp_path, f"({genera})F;", backbone_tips) as report:
        for i in range(n_genera):
            report.taxon_not_monophyletic(f"G{i}", backbone_tips, [f"G{i} a", f"G{i} b"])
            report.placed([f"G{i} b"], "F")
    lines = report.summary().splitlines()
    assert lines[0].startswith(f"Added {n_genera} species; {n_genera} (100%) were placed")
    assert len([line for line in lines if line.startswith("  G")]) == MAX_SUMMARY_TAXA
    assert "  and 2 more taxa" in lines


@pytest.mark.parametrize(
    ("names", "expected"),
    [
        (["a"], "a"),
        (["a", "b"], "a and b"),
        (["a", "b", "c"], "a, b and 1 other"),
        (["a", "b", "c", "d"], "a, b and 2 others"),
    ],
)
def test_fmt_species_list(names, expected):
    assert fmt_species_list(names) == expected
