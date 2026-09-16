import sys

import pytest
from dendropy import Tree

from tact.tree_util import get_tip_labels

execution_number = range(2)


def run_tact(script_runner, tmp_path, config, backbone, args=None):
    if args is None:
        args = []
    config_path = tmp_path / "conf.toml"
    backbone_path = tmp_path / "backbone.tre"
    config_path.write_text(config)
    backbone_path.write_text(backbone)

    result = script_runner.run(
        [
            "tact_add_config",
            "--config",
            str(config_path),
            "--backbone",
            str(backbone_path),
            "--output",
            str(tmp_path / "pytest"),
            "-vv",
            *args,
        ]
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
    new_tips = [f"{focal_clade} tact {x}" for x in range(10)]
    all_tips = {focal_clade, *new_tips}
    mrca_node = res.mrca(taxon_labels=all_tips)
    mrca_tips = {x.taxon.label for x in mrca_node.leaf_iter()}
    assert all_tips == mrca_tips


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("clades", [["A", "B"], ["A", "C"], ["B", "C"]])
def test_poly_singleton(script_runner, execution_number, tmp_path, clades):
    config = f"""
    [[tact]]
    name = "X"
    missing = 10

    [[tact.include]]
    stem = true
    mrca = ["{clades[0]}"]

    [[tact.include]]
    stem = true
    mrca = ["{clades[1]}"]
    """

    backbone = "((A:1,B:1):1,C:2);"
    res = run_tact(script_runner, tmp_path, config, backbone)
    new_tips = set([f"X_tact_{x}" for x in range(10)])

    all_tips = set()
    for clade in clades:
        mrca_node = None
        all_tips.update(clade.mrca)
        tip_node = res.mrca(taxon_labels=all_tips)
        for nn in tip_node.ancestor_iter():
            all_tips.update(get_tip_labels(nn))
            mrca_node = res.mrca(taxon_labels=all_tips)
    raise
