import pytest
import sys
import os

from dendropy import Tree

execution_number = range(2)


def run_tact(script_runner, tmp_path, config, backbone, args=[]):
    config_path = tmp_path / "conf.toml"
    backbone_path = tmp_path / "backbone.tre"
    config_path.write_text(config)
    backbone_path.write_text(backbone)

    result = script_runner.run(
        "tact_add_config",
        "--config",
        config_path,
        "--backbone",
        backbone_path,
        "--output",
        tmp_path / "pytest",
        "-vv",
        *args,
    )
    assert result.returncode == 0

    output = tmp_path / "pytest.newick.tre"
    tacted = Tree.get(path=output, schema="newick", rooting="default-rooted")
    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)
    return tacted


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("focal_clade", ["A", "B", "C"])
def test_lone_singleton(script_runner, execution_number, tmp_path, focal_clade):
    config = f"""
    [[tact]]
    name = "{focal_clade}"
    missing = 10

    [[tact.include]]
    stem = true
    mrca = ["{focal_clade}"]
    """

    backbone = "((A:1,B:1):1,C:2);"
    res = run_tact(script_runner, tmp_path, config, backbone)
    new_tips = [f"{focal_clade}_tact_{x}" for x in range(10)]
    all_tips = set([focal_clade] + new_tips)
    mrca_node = res.mrca(taxon_labels=all_tips)
    mrca_tips = set([x.taxon.label for x in mrca_node.leaf_iter()])
    assert all_tips == mrca_tips
