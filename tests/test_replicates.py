"""Tests for `tact_add_taxa --replicates`."""

import csv
import os
import re

import dendropy
import pytest

from tact import cli_add_taxa, cli_check_trees
from tact.report import PlacementReport
from tact.tree_util import count_newick_trees


@pytest.fixture
def stem2(datadir):
    return [
        "--backbone",
        os.path.join(datadir, "stem2.backbone.tre"),
        "--taxonomy",
        os.path.join(datadir, "stem2.taxonomy.tre"),
    ]


def read(path):
    with open(path, "rb") as fh:
        return fh.read()


def trees(output):
    return dendropy.TreeList.get(path=f"{output}.newick.tre", schema="newick")


def test_replicates(run_cli, stem2, tmp_path):
    output = tmp_path / "out"
    result = run_cli(cli_add_taxa.main, [*stem2, "--output", output, "--replicates", "3", "--seed", "1"])
    forest = trees(output)
    assert len(forest) == 3
    assert all(len(tree.leaf_nodes()) == 25 for tree in forest)
    assert len({tree.as_string(schema="newick") for tree in forest}) == 3
    assert len(dendropy.TreeList.get(path=f"{output}.nexus.tre", schema="nexus")) == 3

    with open(f"{output}.placements.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert {row["replicate"] for row in rows} == {"1", "2", "3"}
    per_replicate = len(rows) // 3
    assert per_replicate > 0 and len(rows) == 3 * per_replicate

    assert result.stdout.startswith(f"Added {per_replicate} species to each of 3 trees")
    # rates are estimated once
    assert read(f"{output}.rates.csv").count(b"\n") > 1


def test_parallel_matches_sequential(run_cli, stem2, tmp_path):
    args = [*stem2, "--replicates", "3", "--seed", "2"]
    run_cli(cli_add_taxa.main, [*args, "--output", tmp_path / "seq"])
    run_cli(cli_add_taxa.main, [*args, "--output", tmp_path / "par", "--cores", "2"])
    for ext in ("newick.tre", "nexus.tre", "placements.csv", "rates.csv"):
        assert read(f"{tmp_path / 'seq'}.{ext}") == read(f"{tmp_path / 'par'}.{ext}")


def test_replicates_depend_on_their_index(run_cli, stem2, tmp_path):
    """Each replicate is seeded from the run's seed and its index, so asking for more trees keeps the first ones."""
    run_cli(cli_add_taxa.main, [*stem2, "--output", tmp_path / "two", "--replicates", "2", "--seed", "3"])
    run_cli(cli_add_taxa.main, [*stem2, "--output", tmp_path / "three", "--replicates", "3", "--seed", "3"])
    two = [tree.as_string(schema="newick") for tree in trees(tmp_path / "two")]
    three = [tree.as_string(schema="newick") for tree in trees(tmp_path / "three")]
    assert three[:2] == two


def test_seed_is_reported_without_seed(run_cli, stem2, tmp_path):
    result = run_cli(cli_add_taxa.main, [*stem2, "--output", tmp_path / "first", "--replicates", "2"])
    match = re.search(r"Made 2 trees; use --seed (\d+) to make the same trees again\.", result.stdout)
    assert match
    again = run_cli(
        cli_add_taxa.main, [*stem2, "--output", tmp_path / "again", "--replicates", "2", "--seed", match.group(1)]
    )
    assert "use --seed" not in again.stdout
    assert read(f"{tmp_path / 'first'}.newick.tre") == read(f"{tmp_path / 'again'}.newick.tre")


def test_one_replicate_keeps_single_tree_output(run_cli, stem2, tmp_path):
    output = tmp_path / "out"
    run_cli(cli_add_taxa.main, [*stem2, "--output", output, "--seed", "4"])
    assert len(trees(output)) == 1
    with open(f"{output}.placements.csv", newline="") as fh:
        assert next(csv.reader(fh))[0] == "species"


def test_report_averages_over_replicates(tmp_path):
    taxonomy = dendropy.Tree.get(data="((G1_a,G1_b,G1_c)G1,(G2_d)G2)F;", schema="newick")
    backbone_tips = {"G1 a", "G2 d"}
    report = PlacementReport(str(tmp_path / "placements.csv"), taxonomy, backbone_tips, replicates=True)
    for replicate in (1, 2):
        one = PlacementReport(None, taxonomy, backbone_tips)
        one.taxon_not_monophyletic("G1", ["G1 a", "G2 d"], ["G1 a", "G1 b", "G1 c"])
        one.placed(["G1 b", "G1 c"] if replicate == 1 else ["G1 b"], "F")
        report.add_replicate(replicate, one.placements)
    report.close()
    lines = report.summary().splitlines()
    assert lines[0].startswith("Added 1.5 species to each of 2 trees; on average 1.5 (100%) were placed")
    assert lines[2] == "  G1: 1.5 species placed in F (1 rank up; broken up by G2 d)"
    with open(tmp_path / "placements.csv", newline="") as fh:
        assert [row[:2] for row in csv.reader(fh)] == [
            ["replicate", "species"],
            ["1", "G1_b"],
            ["1", "G1_c"],
            ["2", "G1_b"],
        ]


@pytest.mark.parametrize("replicates", [1, 2])
def test_check_results_warns_about_replicates(run_cli, stem2, datadir, tmp_path, replicates):
    output = tmp_path / "out"
    run_cli(cli_add_taxa.main, [*stem2, "--output", output, "--replicates", str(replicates), "--seed", "5"])
    args = [
        f"{output}.newick.tre",
        "--backbone",
        os.path.join(datadir, "stem2.backbone.tre"),
        "--taxonomy",
        os.path.join(datadir, "stem2.taxonomy.tre"),
        "--cores",
        "1",
        "--output",
        tmp_path / "check.csv",
    ]
    result = run_cli(cli_check_trees.main, args)
    warning = f"Warning: {output}.newick.tre has 2 trees; checking only the first one."
    assert (warning in result.stderr) == (replicates == 2)


def test_count_newick_trees(tmp_path):
    path = tmp_path / "trees.tre"
    # semicolons in quoted labels and comments do not end a tree
    path.write_text("('a;b',c)[&note=x;y];\n('it''s;',d);\n\n")
    assert count_newick_trees(str(path)) == 2
    path.write_text("(a,b);")
    assert count_newick_trees(str(path)) == 1
