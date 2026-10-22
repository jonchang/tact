"""Tests for passing a taxonomy CSV to `tact_add_taxa` instead of a taxonomy tree."""

import os

import dendropy

from tact import cli_add_taxa, cli_taxonomy

# Binary, ultrametric backbone for the small taxonomies below
BACKBONE = "((A_a:1,A_b:1):1,(B_c:1,B_d:1):1);\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_csv_matches_tree(run_cli, add_taxa, datadir, tmp_path):
    """A CSV gives the same output as the tree `tact_build_taxonomic_tree` builds from it."""
    examples = os.path.join(os.path.dirname(os.path.dirname(datadir)), "examples")
    csv = os.path.join(examples, "Carangaria.csv")
    backbone = os.path.join(examples, "Carangaria.tre")
    tree = tmp_path / "Carangaria.taxonomy.tre"
    run_cli(cli_taxonomy.main, [csv, "--output", tree])

    from_csv = add_taxa(backbone, csv, tmp_path / "csv", seed=0)
    from_tree = add_taxa(backbone, tree, tmp_path / "tree", seed=0)
    assert from_csv == from_tree


def test_csv_problems(run_cli, tmp_path):
    csv = write(tmp_path, "taxonomy.csv", "family,genus,species\nF1,G1,G1 a\nF2,G1,G1 b\nF1,,G1 c\n")
    backbone = write(tmp_path, "backbone.tre", BACKBONE)
    result = run_cli(
        cli_add_taxa.main, ["--taxonomy", csv, "--backbone", backbone, "--output", tmp_path / "out"], exit_code=1
    )
    assert "Found 2 problems in the taxonomy:" in result.output
    assert "* line 3: genus G1 is in F2, but in F1 on line 2" in result.output
    assert "* line 4: empty genus for G1 c (use --fill-empty to fill it in)" in result.output
    assert "Traceback" not in result.output


def test_fill_empty(run_cli, tmp_path):
    csv = write(
        tmp_path,
        "taxonomy.txt",
        "family,genus,species\nF1,A,A a\nF1,A,A b\nF1,A,A e\nF2,B,B c\nF2,B,B d\nF2,,B f\n",
    )
    backbone = write(tmp_path, "backbone.tre", BACKBONE)
    args = ["--taxonomy", csv, "--backbone", backbone, "--output", tmp_path / "out"]
    run_cli(cli_add_taxa.main, args, exit_code=1)

    result = run_cli(cli_add_taxa.main, [*args, "--fill-empty"], seed=0)
    assert "Filled 1 empty cells with placeholder taxa." in result.output
    tree = dendropy.Tree.get(path=str(tmp_path / "out.newick.tre"), schema="newick")
    assert {leaf.taxon.label for leaf in tree.leaf_node_iter()} == {"A a", "A b", "A e", "B c", "B d", "B f"}


def test_fill_empty_needs_csv(run_cli, datadir, tmp_path):
    args = [
        "--taxonomy",
        os.path.join(datadir, "stem.taxonomy.tre"),
        "--backbone",
        os.path.join(datadir, "stem.backbone.tre"),
        "--output",
        tmp_path / "out",
        "--fill-empty",
    ]
    result = run_cli(cli_add_taxa.main, args, exit_code=2)
    assert "--fill-empty only applies to a taxonomy CSV" in result.output


def test_nexus_taxonomy(run_cli, datadir, tmp_path):
    nexus = tmp_path / "taxonomy.nex"
    dendropy.Tree.get(path=os.path.join(datadir, "stem.taxonomy.tre"), schema="newick").write(
        path=str(nexus), schema="nexus"
    )
    args = ["--taxonomy", nexus, "--backbone", os.path.join(datadir, "stem.backbone.tre"), "--output", tmp_path / "o"]
    result = run_cli(cli_add_taxa.main, args, exit_code=2)
    assert "this is a Nexus file" in result.output
