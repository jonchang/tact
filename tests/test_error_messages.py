"""Tests for the errors the command-line tools give for bad input files and options."""

import pytest

from tact import cli_add_taxa, cli_add_toml, cli_check_trees

TAXONOMY = "family,genus,species\nF1,A,A a\nF1,A,A b\nF2,B,B c\nF2,B,B d\n"
BACKBONE = "((A_a:1,A_b:1):1,(B_c:1,B_d:1):1);\n"
NEXUS = "#NEXUS\nbegin trees;\n    tree t = [&R] ((A_a:1,A_b:1):1,(B_c:1,B_d:1):1);\nend;\n"


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


@pytest.fixture
def add_taxa_error(run_cli, tmp_path):
    """Run `tact_add_taxa` on a backbone given as text, expecting it to fail, and return its error output."""

    def run(backbone, *args, output=None, exit_code=2):
        taxonomy = write(tmp_path, "taxonomy.csv", TAXONOMY)
        backbone = write(tmp_path, "backbone.tre", backbone)
        output = output or tmp_path / "out"
        args = ["--taxonomy", taxonomy, "--backbone", backbone, "--output", output, *args]
        result = run_cli(cli_add_taxa.main, args, exit_code=exit_code)
        assert "Traceback" not in result.output
        return result.stderr

    return run


def test_nexus_backbone(add_taxa_error):
    stderr = add_taxa_error(NEXUS)
    assert "Invalid value for '--backbone': this is a Nexus file!" in stderr
    assert "#NEXUS" not in stderr


def test_missing_output_directory(add_taxa_error, tmp_path):
    stderr = add_taxa_error(BACKBONE, output=tmp_path / "nodir" / "out")
    assert f"Invalid value for '--output': directory '{tmp_path / 'nodir'}' does not exist" in stderr


def test_missing_output_directory_checked_first(add_taxa_error, tmp_path):
    """The output directory is checked before the backbone, which is only read once the options are parsed."""
    stderr = add_taxa_error(NEXUS, output=tmp_path / "nodir" / "out")
    assert "'--output'" in stderr
    assert "Nexus" not in stderr


def test_tips_not_in_taxonomy(add_taxa_error):
    tips = ",".join(f"X_{i:02}:1" for i in range(25))
    stderr = add_taxa_error(f"((A_a:1,A_b:1):1,({tips}):1);\n")
    assert "Invalid value for '--backbone': 25 tips are not in the taxonomy:\n* X 00\n* X 01\n" in stderr
    assert "* X 19\n* ... and 5 more\n" in stderr
    assert "X 20" not in stderr
    assert "Add them to the taxonomy, or list outgroups that should stay out of it with --outgroups." in stderr


def test_one_tip_not_in_taxonomy(add_taxa_error):
    stderr = add_taxa_error("((A_a:1,A_b:1):1,(B_c:1,Z_z:1):1);\n")
    assert "1 tip is not in the taxonomy:\n* Z z\nAdd it to the taxonomy" in stderr


def test_outgroups_are_not_missing(run_cli, tmp_path):
    taxonomy = write(tmp_path, "taxonomy.csv", TAXONOMY)
    backbone = write(tmp_path, "backbone.tre", "(((A_a:1,A_b:1):1,(B_c:1,B_d:1):1):1,Z_z:3);\n")
    args = ["--taxonomy", taxonomy, "--backbone", backbone, "--output", tmp_path / "out", "--outgroups", "Z_z"]
    run_cli(cli_add_taxa.main, args)


def test_not_binary(add_taxa_error):
    stderr = add_taxa_error("((A_a:1,A_b:1,B_c:1):1,(B_d:1)X:1);\n")
    assert "the tree is not binary; 2 nodes have other than two branches below them:\n" in stderr
    assert "* 3 branches to A a, A b and 1 other\n* 1 branch to B d\n" in stderr
    assert "at the root" not in stderr


def test_not_binary_at_root(add_taxa_error):
    stderr = add_taxa_error("(A_a:2,A_b:2,(B_c:1,B_d:1):1);\n")
    assert "the tree is not binary; 1 node has other than two branches below it:\n" in stderr
    assert "at the root" in stderr


def test_not_ultrametric(add_taxa_error):
    stderr = add_taxa_error("((A_a:1,A_b:1):1,(B_c:1,B_d:3):1);\n")
    assert (
        "the tree is not ultrametric: tips are between 2.0 (A a) and 4.0 (B d) from the root, "
        "and 3 of 4 tips are closer to it than B d\n"
        "Increase --ultrametricity-precision or use phytools::force.ultrametric in R." in stderr
    )


def test_no_branch_lengths(add_taxa_error):
    stderr = add_taxa_error("((A_a,A_b),(B_c,B_d));\n")
    assert "6 of 6 branches have no length, such as the one to A a; TACT needs a dated tree" in stderr


def test_duplicate_tip(add_taxa_error):
    stderr = add_taxa_error("((A_a:1,A_a:1):1,(B_c:1,B_d:1):1);\n")
    assert "Invalid value for '--backbone': the tip A a appears more than once" in stderr


@pytest.mark.parametrize(
    ("backbone", "message"),
    [
        (" \n", "the file is empty"),
        ("((A_a:1,A_b:1):1,(B_c:1,B_d:1):1)\n", "this is not a Newick tree: it does not end with a semicolon"),
        (
            "((A_a:1,A_b:1):1,(B_c:1,B_d:1):1));\n",
            "this is not a Newick tree: Incomplete or improperly-terminated tree statement "
            "(last character read was ')' instead of a semi-colon ';') (line 1, column 34)",
        ),
    ],
)
def test_not_newick(add_taxa_error, backbone, message):
    assert f"Invalid value for '--backbone': {message}" in add_taxa_error(backbone)


def test_add_config_missing_output_directory(run_cli, tmp_path):
    config = write(
        tmp_path, "config.toml", '[[tact]]\nname = "A"\nmissing = 1\ninclude = [{mrca = ["A_a", "A_b"]}]\n'
    )
    backbone = write(tmp_path, "backbone.tre", BACKBONE)
    args = ["--config", config, "--backbone", backbone, "--output", tmp_path / "nodir" / "out"]
    result = run_cli(cli_add_toml.main, args, exit_code=2)
    assert "Invalid value for '--output': directory" in result.stderr


def test_add_config_nexus_backbone(run_cli, tmp_path):
    config = write(
        tmp_path, "config.toml", '[[tact]]\nname = "A"\nmissing = 1\ninclude = [{mrca = ["A_a", "A_b"]}]\n'
    )
    backbone = write(tmp_path, "backbone.nex", NEXUS)
    result = run_cli(
        cli_add_toml.main, ["--config", config, "--backbone", backbone, "--output", tmp_path / "out"], exit_code=2
    )
    assert "Invalid value for '--backbone': this is a Nexus file" in result.stderr


@pytest.fixture
def check_results_error(run_cli, tmp_path):
    """Run `tact_check_results` with the given files as text, expecting it to fail, and return its error output."""

    def run(simulated=BACKBONE, backbone=BACKBONE, taxonomy="((A_a,A_b)A,(B_c,B_d)B);\n", output=None):
        args = [
            write(tmp_path, "simulated.tre", simulated),
            "--backbone",
            write(tmp_path, "backbone.tre", backbone),
            "--taxonomy",
            write(tmp_path, "taxonomy.tre", taxonomy),
            "--cores",
            1,
        ]
        if output:
            args += ["--output", output]
        result = run_cli(cli_check_trees.main, args, exit_code=2)
        assert "Traceback" not in result.output
        return result.stderr

    return run


def test_check_results_nexus(check_results_error):
    assert "Invalid value for '--backbone': this is a Nexus file" in check_results_error(backbone=NEXUS)
    assert "Invalid value for 'SIMULATED': this is a Nexus file" in check_results_error(simulated=NEXUS)


def test_check_results_taxonomy_csv(check_results_error):
    stderr = check_results_error(taxonomy=TAXONOMY)
    assert "Invalid value for '--taxonomy': this is not a Newick tree, which starts with '('; " in stderr
    assert "make one from a taxonomy CSV with tact_build_taxonomic_tree" in stderr


def test_check_results_missing_output_directory(check_results_error, tmp_path):
    stderr = check_results_error(output=tmp_path / "nodir" / "out.csv")
    assert "Invalid value for '--output'" in stderr
    assert "No such file or directory" in stderr
