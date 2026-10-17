"""Smoke tests that run tact's console scripts in separate processes."""

import os

import pytest
from dendropy import Tree, TreeList

pytestmark = pytest.mark.script_launch_mode("subprocess")

SCRIPTS = ["tact_add_taxa", "tact_add_config", "tact_build_taxonomic_tree", "tact_check_results"]


@pytest.mark.parametrize("script", SCRIPTS)
def test_help(script_runner, script):
    result = script_runner.run([script, "--help"])
    assert result.returncode == 0
    assert result.stdout.startswith(f"Usage: {script} ")


def test_build_add_and_check(script_runner, tmp_path):
    taxonomy_csv = tmp_path / "taxonomy.csv"
    taxonomy_csv.write_text(
        "family,genus,species\n"
        "F1,G1,G1 a\nF1,G1,G1 b\nF1,G1,G1 c\nF1,G2,G2 d\nF1,G2,G2 e\n"
        "F2,G3,G3 f\nF2,G3,G3 g\nF2,G4,G4 h\n"
    )
    backbone = tmp_path / "backbone.tre"
    backbone.write_text("(((G1_a:3,G2_d:3):2,G2_e:5):5,(G3_f:4,G4_h:4):6);")
    taxonomy = tmp_path / "taxonomy.tre"
    output = tmp_path / "out"

    result = script_runner.run(["tact_build_taxonomic_tree", str(taxonomy_csv), "--output", str(taxonomy)])
    assert result.returncode == 0

    result = script_runner.run(
        ["tact_add_taxa", "--backbone", str(backbone), "--taxonomy", str(taxonomy), "--output", str(output)]
    )
    assert result.returncode == 0
    tacted = Tree.get(path=f"{output}.newick.tre", schema="newick")
    assert len(tacted.leaf_nodes()) == 8
    assert os.path.exists(f"{output}.rates.csv")

    result = script_runner.run(
        [
            "tact_check_results",
            f"{output}.newick.tre",
            "--backbone",
            str(backbone),
            "--taxonomy",
            str(taxonomy),
            "--output",
            str(tmp_path / "check.csv"),
            "--cores=1",
        ]
    )
    assert result.returncode == 0


def test_add_config(script_runner, tmp_path):
    config = tmp_path / "config.toml"
    config.write_text('[[tact]]\nname = "A"\nmissing = 3\n\n[[tact.include]]\nmrca = ["A", "B"]\n')
    backbone = tmp_path / "backbone.tre"
    backbone.write_text("((A:1,B:1):1,C:2);")
    output = tmp_path / "out"

    result = script_runner.run(
        [
            "tact_add_config",
            "--config",
            str(config),
            "--backbone",
            str(backbone),
            "--output",
            str(output),
            "--replicates=2",
            "--cores=1",
        ]
    )
    assert result.returncode == 0
    trees = TreeList.get(path=f"{output}.newick.tre", schema="newick")
    assert len(trees) == 2
    assert all(len(tree.leaf_nodes()) == 6 for tree in trees)
