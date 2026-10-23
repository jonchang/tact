"""Tests for names that differ only in spaces, underscores or case, which TACT treats as the same name."""

import csv

import pytest

from tact import cli_add_taxa, cli_add_toml, cli_check_trees
from tact.cli_taxonomy import TaxonomyError, read_taxonomy
from tact.names import clean_name, name_key, written_name

# A c is missing from the backbone
TAXONOMY = """family,genus,species
F1,A,A a
F1,A,A b
F1,A,A c
F1,B,B c
F1,B,B d
"""
BACKBONE = "((A_a:1,A_b:1):1,(B_c:1,B_d:1):1);\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def add_taxa(run_cli, tmp_path, taxonomy, backbone, *args, exit_code=0):
    """Run `tact_add_taxa` on a taxonomy and backbone given as text."""
    taxonomy_path = write(tmp_path, "taxonomy.csv" if not taxonomy.startswith("(") else "taxonomy.tre", taxonomy)
    backbone_path = write(tmp_path, "backbone.tre", backbone)
    args = ["--taxonomy", taxonomy_path, "--backbone", backbone_path, "--output", tmp_path / "out", *args]
    result = run_cli(cli_add_taxa.main, args, seed=0, exit_code=exit_code)
    assert "Traceback" not in result.output
    return result


@pytest.mark.parametrize(
    ("name", "clean"),
    [
        ("Gobius_niger", "Gobius niger"),
        ("Gobius__niger ", "Gobius niger"),
        ("Gobius" + chr(0xA0) + "niger", "Gobius niger"),
    ],
)
def test_clean_name(name, clean):
    assert clean_name(name) == clean
    assert name_key(name) == name_key("gobius NIGER")
    assert written_name(clean) == "Gobius_niger"


def test_hyphens_are_kept():
    assert clean_name("Quietula_y-cauda") == "Quietula y-cauda"
    assert written_name("Quietula y-cauda") == "Quietula_y-cauda"
    assert name_key("Quietula y-cauda") != name_key("Quietula ycauda")
    assert name_key("Quietula y-cauda") != name_key("Quietula y cauda")


def test_quietula_y_cauda(run_cli, tmp_path):
    """Hyphenated species names match between a taxonomy CSV and a backbone, and are written unquoted."""
    taxonomy = (
        "family,genus,species\n"
        "Gobiidae,Quietula,Quietula y-cauda\n"
        "Gobiidae,Quietula,Quietula guaymasiae\n"
        "Gobiidae,Gobius,Gobius niger\n"
        "Gobiidae,Gobius,Gobius paganellus\n"
        "Gobiidae,Gobius,Gobius x-ray\n"
    )
    backbone = "((Quietula_y-cauda:1,Quietula_guaymasiae:1):1,(Gobius_niger:1,Gobius_paganellus:1):1);\n"
    add_taxa(run_cli, tmp_path, taxonomy, backbone)
    tree = (tmp_path / "out.newick.tre").read_text()
    assert "Quietula_y-cauda:" in tree
    assert "Gobius_x-ray:" in tree
    assert "'" not in tree
    assert [row["species"] for row in read_csv(tmp_path / "out.placements.csv")] == ["Gobius_x-ray"]


def test_csv_with_underscores_matches_backbone(run_cli, tmp_path):
    add_taxa(run_cli, tmp_path, TAXONOMY.replace(" ", "_"), BACKBONE)
    # names are written with underscores, unquoted, in every file
    tree = (tmp_path / "out.newick.tre").read_text()
    assert "A_c:" in tree
    assert "'" not in tree
    (row,) = read_csv(tmp_path / "out.placements.csv")
    assert row["species"] == "A_c"
    assert {row["taxon"] for row in read_csv(tmp_path / "out.rates.csv")} == {"__TAXONOMIC_ROOT__", "F1", "A", "B"}


def test_newick_taxonomy_with_quoted_underscores(run_cli, tmp_path):
    """A quoted label keeps its underscores in DendroPy, but TACT still reads them as spaces."""
    taxonomy = "((('A_a','A_b','A_c')A,('B_c','B_d')B)F1);\n"
    result = add_taxa(run_cli, tmp_path, taxonomy, BACKBONE)
    assert "Added 1 species" in result.stdout


def test_csv_species_spelled_two_ways(tmp_path):
    path = write(tmp_path, "taxonomy.csv", TAXONOMY + "F1,A,A_a\nF1,A,a a\nF1,A,a a\n")
    with pytest.raises(TaxonomyError) as info:
        read_taxonomy(str(path))
    assert info.value.problems == [
        "line 7: A_a is spelled A a on line 2; spell each species the same way on every line",
        "line 8: a a is spelled A a on line 2; spell each species the same way on every line",
    ]


def test_csv_rank_spelled_two_ways(run_cli, tmp_path):
    taxonomy = "family,genus,species\nF 1,A,A a\nF_1,A,A b\nf__1,A,A c\nF 1,B,B c\nF_1,B,B d\n"
    parsed = read_taxonomy(str(write(tmp_path, "t.csv", taxonomy)))
    assert list(parsed.root.children) == ["F 1"]
    assert parsed.respelled == {"F 1": [("F 1", 2), ("F_1", 3), ("f__1", 4)]}
    result = add_taxa(run_cli, tmp_path, taxonomy, BACKBONE)
    assert (
        "Note: these names are spelled more than one way, and each is read as the same taxon:\n"
        "F 1: F 1 (line 2), F_1 (line 3), f__1 (line 4)\n" in result.stdout
    )


def test_backbone_tips_spelled_two_ways(run_cli, tmp_path):
    backbone = "((A_a:1,'A a':1):1,(B_c:1,b__C:1):1);\n"
    result = add_taxa(run_cli, tmp_path, TAXONOMY, backbone, exit_code=2)
    assert (
        "Invalid value for '--backbone': 2 names are spelled more than one way, "
        "differing only in spaces, underscores or case:\n"
        "* tips A_a and A a\n"
        "* tips B_c and b__C\n" in result.stderr
    )


def test_taxonomy_tree_taxa_spelled_two_ways(run_cli, tmp_path):
    taxonomy = "(((A_a,A_b,A_c)Genus_A,(B_c,B_d)'genus A')F1);\n"
    result = add_taxa(run_cli, tmp_path, taxonomy, BACKBONE, exit_code=2)
    assert "Invalid value for '--taxonomy': 1 name is spelled" in result.stderr
    assert "* taxa Genus_A and genus A\n" in result.stderr


def test_backbone_case_differs_from_taxonomy(run_cli, tmp_path):
    result = add_taxa(run_cli, tmp_path, TAXONOMY, "((a_A:1,A_b:1):1,(B_c:1,B_d:1):1);\n")
    assert "Note: these backbone tips are spelled differently in the taxonomy, which TACT uses:\n* a A => A a\n" in (
        result.stdout
    )
    assert "A_a:" in (tmp_path / "out.newick.tre").read_text()


def test_outgroups_with_underscores_or_spaces(run_cli, tmp_path):
    backbone = "((((A_a:1,A_b:1):1,(B_c:1,B_d:1):1):1,Out_one:3):1,Out_two:4);\n"
    add_taxa(run_cli, tmp_path, TAXONOMY, backbone, "--outgroups", "Out one,Out__two")
    assert "Out_one:" in (tmp_path / "out.newick.tre").read_text()


def test_add_config_mrca_with_spaces(run_cli, tmp_path):
    config = write(
        tmp_path, "config.toml", '[[tact]]\nname = "A_new"\nmissing = 1\ninclude = [{mrca = ["A a", "A_b"]}]\n'
    )
    backbone = write(tmp_path, "backbone.tre", BACKBONE)
    args = ["--config", config, "--backbone", backbone, "--output", tmp_path / "out", "--replicates", 1, "--cores", 1]
    run_cli(cli_add_toml.main, args, seed=0)
    assert "A_new_tact_0:" in (tmp_path / "out.newick.tre").read_text()


def test_check_results_matches_names(run_cli, tmp_path):
    """Trees that spell names differently are compared by name, and the report spells taxa with underscores."""
    simulated = write(tmp_path, "simulated.tre", "(((A_a:1,A_b:1):0.5,'A c':1.5):0.5,(B_c:1,B_d:1):1);\n")
    backbone = write(tmp_path, "backbone.tre", "(('A a':1,A_b:1):1,(B_c:1,B_d:1):1);\n")
    taxonomy = write(tmp_path, "taxonomy.tre", "(((A_a,A_b,A_c)Genus_A,(B_c,B_d)Genus_B)F1);\n")
    output = tmp_path / "check.csv"
    args = [simulated, "--backbone", backbone, "--taxonomy", taxonomy, "--output", output, "--cores", 1]
    run_cli(cli_check_trees.main, args)
    rows = {row["node"]: row for row in read_csv(output)}
    assert rows.keys() == {"F1", "Genus_A", "Genus_B"}
    assert rows["Genus_A"]["taxonomy_tips"] == "3"
    assert rows["Genus_A"]["backbone_tips"] == "2"
    assert rows["Genus_A"]["simulated_tips"] == "3"
