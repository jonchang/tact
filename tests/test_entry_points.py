"""Smoke tests that run tact's console scripts in separate processes."""

import os
import subprocess
import sys

import pytest
from dendropy import Tree, TreeList

pytestmark = pytest.mark.script_launch_mode("subprocess")

SCRIPTS = ["tact_add_taxa", "tact_add_config", "tact_build_taxonomic_tree", "tact_check_results"]
CLI_NUMPY_MODULES = ["tact.cli_add_taxa", "tact.cli_add_toml", "tact.cli_check_trees"]
BLAS_THREAD_VARS = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")

# Imports a module and prints the BLAS thread variables as they were when NumPy started loading
BLAS_THREADS_AT_NUMPY_IMPORT = """
import importlib, os, sys

class Spy:
    def find_spec(self, name, path=None, target=None):
        if name == "numpy":
            sys.meta_path.remove(self)
            print(",".join(os.environ.get(var, "") for var in sys.argv[2:]))

sys.meta_path.insert(0, Spy())
importlib.import_module(sys.argv[1])
"""


def blas_threads_at_numpy_import(module, env):
    result = subprocess.run(
        [sys.executable, "-c", BLAS_THREADS_AT_NUMPY_IMPORT, module, *BLAS_THREAD_VARS],
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip().split(",")


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


@pytest.mark.parametrize("module", CLI_NUMPY_MODULES)
def test_blas_threads_limited_before_numpy_loads(module):
    env = {key: value for key, value in os.environ.items() if key not in BLAS_THREAD_VARS}
    assert blas_threads_at_numpy_import(module, env) == ["1", "1", "1"]


@pytest.mark.parametrize("module", CLI_NUMPY_MODULES)
def test_blas_threads_keep_user_setting(module):
    env = {**os.environ, "OPENBLAS_NUM_THREADS": "4", "OMP_NUM_THREADS": "2", "MKL_NUM_THREADS": "3"}
    assert blas_threads_at_numpy_import(module, env) == ["4", "2", "3"]
