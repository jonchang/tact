import pytest
import sys
import os

from dendropy import Tree

execution_number = range(5)

def run_tact(script_runner, datadir, stem):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    taxed = Tree.get(path=taxonomy, schema="newick")
    bbone = Tree.get(path=backbone, schema="newick")
    result = script_runner.run("tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", stem, "-vv")
    assert result.returncode == 0
    output = stem + ".newick.tre"
    tacted = Tree.get(path=output, schema="newick")
    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)
    result = script_runner.run("tact_check_results", output, "--taxonomy", taxonomy, "--backbone", backbone, "--output", stem + ".check.csv", "--cores=1")
    assert result.returncode == 0
    return (tacted, taxed, bbone)

@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "intrusion", "short_branch", "stem"])
@pytest.mark.script_launch_mode('subprocess')
def test_monophyly(script_runner, execution_number, datadir, stem):
    tacted, taxed, bbone = run_tact(script_runner, datadir, stem)
    extant = set([x.taxon.label for x in bbone.leaf_nodes()])
    for node in taxed.postorder_internal_node_iter(exclude_seed_node=True):
        label = node.label
        expected = set([x.taxon.label for x in node.leaf_nodes()])
        our_extant = extant & expected
        if len(our_extant) > 0:
            bbone_node = bbone.mrca(taxon_labels=our_extant)
            bbone_tips = set([x.taxon.label for x in bbone_node.leaf_nodes()])
            if bbone_tips != our_extant:
                continue
        mrca = tacted.mrca(taxon_labels=expected)
        actual = set([x.taxon.label for x in mrca.leaf_nodes()])
        assert expected == actual

@pytest.mark.parametrize('execution_number', execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "intrusion", "short_branch"])
@pytest.mark.script_launch_mode('subprocess')
def test_short_branch(script_runner, execution_number, datadir, stem):
    tacted, taxed, bbone = run_tact(script_runner, datadir, stem)
    n_short = 0
    for leaf in tacted.leaf_node_iter():
        if leaf.edge.length < 0.1:
            n_short += 1
    assert n_short <= 15


@pytest.mark.parametrize('execution_number', execution_number)
@pytest.mark.parametrize("stem", ["weirdness"])
@pytest.mark.script_launch_mode('subprocess')
def test_stem_clade_attachment(script_runner, execution_number, datadir, stem):
    tacted, taxed, bbone = run_tact(script_runner, datadir, "weirdness")
