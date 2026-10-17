import os
import random

import numpy as np
import pytest
from click.testing import CliRunner

from tact import cli_add_taxa, console
from tact.cli_taxonomy import build_taxonomic_tree


@pytest.fixture
def datadir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@pytest.fixture
def run_cli():
    """Run a TACT command in-process and return a `click.testing.Result`.

    Invoke it as `run_cli(command, args, seed=None, exit_code=0)`.

    The commands keep state at module level, so each call seeds the random number generators when given a
    `seed` and performs necessary cleanups such as emptying shared caches and removing log handlers.
    """

    def run(command, args, seed=None, exit_code=0):
        if seed is not None:
            random.seed(seed)
            # SciPy's annealing fallback draws from NumPy's global state
            np.random.seed(seed)  # noqa: NPY002
        cli_add_taxa.mrca_rates.clear()
        try:
            result = CliRunner().invoke(command, [str(arg) for arg in args], catch_exceptions=False)
        finally:
            console.reset_logging()
        assert result.exit_code == exit_code, result.output
        return result

    return run


@pytest.fixture
def add_taxa(run_cli):
    """Run `tact_add_taxa` in-process and return its tree and rate files as bytes, keyed by extension."""

    def run(backbone, taxonomy, output, *args, seed=None):
        run_cli(cli_add_taxa.main, ["--backbone", backbone, "--taxonomy", taxonomy, "--output", output, *args], seed)
        outputs = {}
        for ext in ("newick.tre", "nexus.tre", "rates.csv"):
            with open(f"{output}.{ext}", "rb") as fh:
                outputs[ext] = fh.read()
        return outputs

    return run


@pytest.fixture
def carangaria(datadir, tmp_path):
    """Paths to the Carangaria example backbone and a taxonomy tree built from its CSV."""
    examples = os.path.join(os.path.dirname(os.path.dirname(datadir)), "examples")
    taxonomy = tmp_path / "Carangaria.taxonomy.tre"
    build_taxonomic_tree(os.path.join(examples, "Carangaria.csv")).write_to_path(str(taxonomy), schema="newick")
    return os.path.join(examples, "Carangaria.tre"), taxonomy


@pytest.fixture
def ages():
    return [
        20.934955,
        17.506532,
        16.64467,
        15.380987,
        14.547092000000001,
        13.664577999999999,
        13.289480000000001,
        11.667099,
        9.799231,
        9.510413,
        9.029556000000001,
        8.806255,
        8.770727,
        8.480102,
        6.984475999999999,
        6.706684,
        2.11319,
        0.545689,
        0.147482,
    ]


@pytest.fixture
def sampling():
    return 0.869565217391


@pytest.fixture
def birth():
    return 0.09115578881894915


@pytest.fixture
def death():
    return 0.0
