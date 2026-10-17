"""Tests for log messages on the console and in the log file."""

import logging
import os

import pytest

from tact import cli_add_taxa, cli_add_toml, console

INFO_MESSAGE = "Backbone needs to add"
DEBUG_MESSAGE = "Computing root birth and death rates."


@pytest.fixture
def stem2(datadir):
    return os.path.join(datadir, "stem2.backbone.tre"), os.path.join(datadir, "stem2.taxonomy.tre")


def add_taxa_args(stem2, output, *args):
    backbone, taxonomy = stem2
    return ["--backbone", backbone, "--taxonomy", taxonomy, "--output", output, *args]


def read_log(output):
    with open(f"{output}.log.txt", encoding="utf-8") as fh:
        return fh.read()


@pytest.mark.parametrize(
    ("verbose", "console_info", "console_debug", "file_debug"),
    [
        ([], False, False, False),
        (["-v"], True, False, False),
        (["-vv"], True, True, True),
    ],
)
def test_log_levels(run_cli, stem2, tmp_path, verbose, console_info, console_debug, file_debug):
    output = tmp_path / "out"
    result = run_cli(cli_add_taxa.main, add_taxa_args(stem2, output, *verbose))
    assert (INFO_MESSAGE in result.stderr) == console_info
    assert (DEBUG_MESSAGE in result.stderr) == console_debug
    log = read_log(output)
    # the log file always has INFO messages
    assert INFO_MESSAGE in log
    assert (DEBUG_MESSAGE in log) == file_debug


def test_add_config_logs_to_console(run_cli, tmp_path):
    config = tmp_path / "config.toml"
    config.write_text('[[tact]]\nname = "A"\nmissing = 3\n\n[[tact.include]]\nmrca = ["A", "B"]\n')
    backbone = tmp_path / "backbone.tre"
    backbone.write_text("((A:1,B:1):1,C:2);")
    output = tmp_path / "out"
    args = ["--config", config, "--backbone", backbone, "--output", output, "--replicates=1", "--cores=1"]
    result = run_cli(cli_add_toml.main, [*args, "-v"])
    assert "Reading configuration" in result.stderr
    assert "Reading configuration" in read_log(output)


def test_setup_logging_replaces_handlers(tmp_path):
    log = logging.getLogger("tact.test_console")
    try:
        console.setup_logging(str(tmp_path / "first"), 0)
        console.setup_logging(str(tmp_path / "second"), 0)
        log.info("only in the second log")
    finally:
        console.reset_logging()
    assert read_log(tmp_path / "first") == ""
    assert read_log(tmp_path / "second") == "only in the second log\n"
