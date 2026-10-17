import os
import sys

import pytest
from dendropy import Tree

from tact import cli_add_taxa, cli_check_trees

execution_number = range(2)


def run_tact(run_cli, datadir, stem, tmp_path, *args):
    """Add taxa to a test backbone, check the result with `tact_check_results`, and return the trees."""
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    taxed = Tree.get(path=taxonomy, schema="newick")
    bbone = Tree.get(path=backbone, schema="newick", rooting="default-rooted")
    output_base = tmp_path / f"tact-pytest-{stem}"
    add_args = ["--taxonomy", taxonomy, "--backbone", backbone, "--output", output_base]
    run_cli(cli_add_taxa.main, [*add_args, "-vv", *args])
    output = str(output_base) + ".newick.tre"
    tacted = Tree.get(path=output, schema="newick", rooting="default-rooted")
    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)
    check_args = [output, "--taxonomy", taxonomy, "--backbone", backbone, "--output", f"{output_base}.check.csv"]
    run_cli(cli_check_trees.main, [*check_args, "--cores=1"])
    return (tacted, taxed, bbone)


@pytest.mark.parametrize("execution_number", execution_number)
def test_yule(run_cli, execution_number, datadir, tmp_path):
    run_tact(run_cli, datadir, "stem2", tmp_path, "--yule")


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "intrusion", "short_branch", "stem"])
def test_monophyly(run_cli, execution_number, datadir, stem, tmp_path):
    tacted, taxed, bbone = run_tact(run_cli, datadir, stem, tmp_path)
    extant = {x.taxon.label for x in bbone.leaf_nodes()}
    for node in taxed.postorder_internal_node_iter(exclude_seed_node=True):
        expected = {x.taxon.label for x in node.leaf_nodes()}
        our_extant = extant & expected
        if len(our_extant) > 0:
            bbone_node = bbone.mrca(taxon_labels=our_extant)
            bbone_tips = {x.taxon.label for x in bbone_node.leaf_nodes()}
            if bbone_tips != our_extant:
                continue
        mrca = tacted.mrca(taxon_labels=expected)
        actual = {x.taxon.label for x in mrca.leaf_nodes()}
        assert expected == actual


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "short_branch"])
def test_short_branch(run_cli, execution_number, datadir, stem, tmp_path):
    tacted, _taxed, _bbone = run_tact(run_cli, datadir, stem, tmp_path)
    n_short = 0
    for leaf in tacted.leaf_node_iter():
        if leaf.edge.length < 0.1:
            n_short += 1
    assert n_short <= 20


@pytest.mark.parametrize("execution_number", execution_number)
def test_stem_clade_attachment(run_cli, execution_number, datadir, tmp_path):
    tacted, _taxed, _bbone = run_tact(run_cli, datadir, "stem2", tmp_path)
    tacted.calc_node_ages()
    tacted.update_bipartitions()
    node = tacted.mrca(taxon_labels=["c1", "c2", "c3", "c4", "c5"])
    ages = [x.age < 15.16 for x in node.postorder_internal_node_iter()]
    assert any(ages)
