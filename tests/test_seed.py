"""Check that `--seed` makes runs reproducible across processes with different hash seeds."""

import os

import pytest

pytestmark = pytest.mark.script_launch_mode("subprocess")

EXAMPLES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples")


def run_seeded(script_runner, args, output, seed, hash_seed):
    """Run a command with `--seed` under `PYTHONHASHSEED` and return its output files as bytes, keyed by extension."""
    env = {**os.environ, "PYTHONHASHSEED": str(hash_seed)}
    result = script_runner.run([*args, "--output", str(output), "--seed", str(seed)], env=env)
    assert result.returncode == 0, result.stderr
    outputs = {}
    for path in output.parent.iterdir():
        if path.name.startswith(output.name + ".") and not path.name.endswith(".log"):
            outputs[path.name.removeprefix(output.name + ".")] = path.read_bytes()
    return outputs


def test_add_taxa_seed(script_runner, tmp_path):
    args = [
        "tact_add_taxa",
        "--backbone",
        os.path.join(EXAMPLES, "Carangaria.tre"),
        "--taxonomy",
        os.path.join(EXAMPLES, "Carangaria.csv"),
    ]
    first = run_seeded(script_runner, args, tmp_path / "first", seed=42, hash_seed=1)
    again = run_seeded(script_runner, args, tmp_path / "again", seed=42, hash_seed=2)
    other = run_seeded(script_runner, args, tmp_path / "other", seed=43, hash_seed=1)

    for ext in ("newick.tre", "rates.csv", "placements.csv"):
        assert first[ext] == again[ext], ext
    assert first["newick.tre"] != other["newick.tre"]


def test_add_config_seed(script_runner, tmp_path):
    config = tmp_path / "config.toml"
    config.write_text(
        '[[tact]]\nname = "AB"\nmissing = 5\n\n[[tact.include]]\nmrca = ["A", "B"]\n\n'
        '[[tact]]\nname = "CD"\nmissing = 5\n\n[[tact.include]]\nmrca = ["C", "D"]\n'
    )
    backbone = tmp_path / "backbone.tre"
    backbone.write_text("(((A:1,B:1):2,(C:2,D:2):1):1,E:4);")
    args = ["tact_add_config", "--config", str(config), "--backbone", str(backbone), "--replicates=2", "--cores=2"]
    first = run_seeded(script_runner, args, tmp_path / "first", seed=42, hash_seed=1)
    again = run_seeded(script_runner, args, tmp_path / "again", seed=42, hash_seed=2)
    other = run_seeded(script_runner, args, tmp_path / "other", seed=43, hash_seed=1)

    assert first["newick.tre"] == again["newick.tre"]
    assert first["newick.tre"] != other["newick.tre"]
    # each replicate gets its own seed
    replicates = first["newick.tre"].decode().splitlines()
    assert len(replicates) == 2
    assert replicates[0] != replicates[1]
