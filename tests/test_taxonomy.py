"""Tests for building taxonomic trees from CSV files."""

import random
import time

import dendropy
import pytest

from tact import cli_taxonomy
from tact.cli_taxonomy import ROOT, TaxonomyError, build_taxonomic_tree, quote, read_taxonomy


def write_csv(tmp_path, text, name="taxonomy.csv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


def clades(tree):
    """Map each labelled internal node to the tip labels below it."""
    out = {}
    for node in tree.postorder_internal_node_iter():
        if node.label:
            out[node.label] = {leaf.taxon.label for leaf in node.leaf_iter()}
    return out


def problems(tmp_path, text, **kwargs):
    with pytest.raises(TaxonomyError) as info:
        read_taxonomy(write_csv(tmp_path, text), **kwargs)
    return info.value.problems


def test_builds_tree(tmp_path):
    tree = build_taxonomic_tree(
        write_csv(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G1,G1 b\nF1,G2,G2 c\nF2,G3,G3 d\n")
    )
    assert clades(tree) == {
        ROOT: {"G1 a", "G1 b", "G2 c", "G3 d"},
        "F1": {"G1 a", "G1 b", "G2 c"},
        "F2": {"G3 d"},
        "G1": {"G1 a", "G1 b"},
        "G2": {"G2 c"},
        "G3": {"G3 d"},
    }


def test_row_order_does_not_matter(tmp_path):
    rows = ["F1,G1,G1 a", "F1,G1,G1 b", "F1,G2,G2 c", "F2,G3,G3 d", "F2,G3,G3 e"]
    shuffled = rows[:]
    random.Random(1).shuffle(shuffled)
    first = read_taxonomy(write_csv(tmp_path, "family,genus,species\n" + "\n".join(rows) + "\n", "a.csv"))
    second = read_taxonomy(write_csv(tmp_path, "family,genus,species\n" + "\n".join(shuffled) + "\n", "b.csv"))
    assert first.to_newick() == second.to_newick()


def test_genus_in_two_families(tmp_path):
    found = problems(tmp_path, "family,genus,species\nF1,G1,G1 a\nF2,G1,G1 b\nF2,G1,G1 c\n")
    # reported once, not once per row
    assert found == ["line 3: genus G1 is in F2, but in F1 on line 2"]


def test_wrong_number_of_cells(tmp_path):
    found = problems(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G1\nF2,G2,G2 b,extra\n")
    assert found == ["line 3: expected 3 cells, found 2", "line 4: expected 3 cells, found 4"]


def test_empty_cells(tmp_path):
    text = "order,family,genus,species\nO1,F1,G1,G1 a\nO1,,G2,G2 b\n,,G3,G3 c\nO1,F1,G1,\n"
    found = problems(tmp_path, text)
    assert found == [
        "line 3: empty family for G2 b (use --fill-empty to fill it in)",
        "line 4: empty order, family for G3 c (use --fill-empty to fill it in)",
        "line 5: no species name",
    ]


def test_fill_empty(tmp_path):
    text = "order,family,genus,species\nO1,F1,G1,G1 a\nO1,,G2,G2 b\nO1,,G2,G2 c\n,,G3,G3 d\n"
    taxonomy = read_taxonomy(write_csv(tmp_path, text), fill_empty=True)
    assert taxonomy.filled == 4
    tree = dendropy.Tree.get(data=taxonomy.to_newick(), schema="newick")
    found = clades(tree)
    # placeholders are named after the nearest named taxon below them
    assert found["G2 family"] == {"G2 b", "G2 c"}
    assert found["G3 order"] == found["G3 family"] == {"G3 d"}
    assert found["O1"] == {"G1 a", "G2 b", "G2 c"}


def test_name_at_two_ranks(tmp_path):
    # Cichla is both a genus and its nominal subgenus; the lower rank is renamed whatever the row order
    for order in ([0, 1], [1, 0]):
        rows = ["Cichlidae,Cichla,Cichla,Cichla temensis", "Cichlidae,Cichla,Other,Cichla other"]
        text = "family,genus,subgenus,species\n" + "\n".join(rows[i] for i in order) + "\n"
        taxonomy = read_taxonomy(write_csv(tmp_path, text))
        assert taxonomy.renamed == {("Cichla", "subgenus"): "Cichla subgenus"}
        found = clades(dendropy.Tree.get(data=taxonomy.to_newick(), schema="newick"))
        assert found["Cichla subgenus"] == {"Cichla temensis"}
        assert found["Cichla"] == {"Cichla temensis", "Cichla other"}


def test_species_names_are_not_renamed(tmp_path):
    taxonomy = read_taxonomy(write_csv(tmp_path, "family,genus,species\nF1,Foo,Foo\n"))
    assert taxonomy.renamed == {("Foo", "genus"): "Foo genus"}
    tree = dendropy.Tree.get(data=taxonomy.to_newick(), schema="newick")
    assert [leaf.taxon.label for leaf in tree.leaf_node_iter()] == ["Foo"]


def test_duplicate_species(tmp_path):
    # an identical row is ignored
    taxonomy = read_taxonomy(write_csv(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G1,G1 a\n", "same.csv"))
    assert list(taxonomy.root.children["F1"].children["G1"].children) == ["G1 a"]
    found = problems(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G2,G1 a\n")
    assert found == ["line 3: G1 a is also listed on line 2, under different ranks"]


def test_header_problems(tmp_path):
    assert problems(tmp_path, "genus,genus,species\nG1,G1,G1 a\n") == ["line 1: duplicated column names: genus"]
    assert problems(tmp_path, "species\nG1 a\n") == [
        "line 1: there must be at least one rank column and a species column"
    ]
    assert problems(tmp_path, "") == ["line 1: the file is empty"]


def test_byte_order_mark_and_whitespace(tmp_path):
    taxonomy = read_taxonomy(write_csv(tmp_path, "﻿family , genus,species\n F1,G1 , G1 a \n\n"))
    assert taxonomy.ranks == ["family", "genus"]
    assert list(taxonomy.root.children["F1"].children["G1"].children) == ["G1 a"]


@pytest.mark.parametrize(
    "label",
    ["Gobius niger", "Quietula y-cauda", "Gobius (Gobius) niger", "Gobius niger 'A'", "Gobius sp. A:B", 'say "hi"'],
)
def test_labels_survive_newick(tmp_path, label):
    text = f'family,genus,species\nF1,G1,"{label.replace(chr(34), chr(34) * 2)}"\n'
    tree = build_taxonomic_tree(write_csv(tmp_path, text))
    assert [leaf.taxon.label for leaf in tree.leaf_node_iter()] == [label]
    assert quote("plain name") == "plain_name"


def test_many_problems_are_summarized(tmp_path):
    rows = "".join(f"F1,,G{i} a\n" for i in range(30))
    with pytest.raises(TaxonomyError) as info:
        read_taxonomy(write_csv(tmp_path, "family,genus,species\n" + rows))
    message = info.value.message
    assert message.startswith("Found 30 problems in the taxonomy:")
    assert message.endswith("* ... and 10 more")


def test_large_taxonomy_is_fast(tmp_path):
    # 50,000 species in 5,000 genera; the previous builder searched the whole tree for each new genus
    rows = [f"O{g // 1000},F{g // 20},G{g},G{g} s{s}" for g in range(5000) for s in range(10)]
    random.Random(0).shuffle(rows)
    path = write_csv(tmp_path, "order,family,genus,species\n" + "\n".join(rows) + "\n")
    start = time.perf_counter()
    newick = read_taxonomy(path).to_newick()
    assert time.perf_counter() - start < 5
    assert newick.count(",") == 50000 - 1


def test_cli(run_cli, tmp_path):
    path = write_csv(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G1,G1 b\n")
    output = tmp_path / "out.tre"
    result = run_cli(cli_taxonomy.main, [path, "--output", output])
    assert result.stdout == f"Output written to: {output}\n"
    tree = dendropy.Tree.get(path=str(output), schema="newick")
    assert clades(tree)["G1"] == {"G1 a", "G1 b"}


def test_cli_nexus(run_cli, tmp_path):
    path = write_csv(tmp_path, "family,genus,species\nF1,G1,G1 a\nF1,G1,G1 b\n")
    output = tmp_path / "out.nex"
    run_cli(cli_taxonomy.main, [path, "--output", output, "--schema", "nexus"])
    tree = dendropy.Tree.get(path=str(output), schema="nexus")
    assert clades(tree)["G1"] == {"G1 a", "G1 b"}


def test_cli_reports_problems(run_cli, tmp_path):
    path = write_csv(tmp_path, "family,genus,species\nF1,G1,G1 a\nF2,G1,G1 b\nF1,,G2 c\n")
    result = run_cli(cli_taxonomy.main, [path, "--output", tmp_path / "out.tre"], exit_code=1)
    assert "Found 2 problems in the taxonomy:" in result.stderr
    assert "line 3: genus G1 is in F2, but in F1 on line 2" in result.stderr
    assert "line 4: empty genus for G2 c" in result.stderr
    assert not (tmp_path / "out.tre").exists()
