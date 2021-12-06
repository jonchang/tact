import pytest
import sys
import os

from dendropy import Tree

execution_number = range(2)


def run_tact(script_runner, datadir, stem):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    taxed = Tree.get(path=taxonomy, schema="newick")
    bbone = Tree.get(path=backbone, schema="newick", rooting="default-rooted")
    result = script_runner.run(
        "tact_add_taxa", "--taxonomy", taxonomy, "--backbone", backbone, "--output", ".tact-pytest-" + stem, "-vv"
    )
    assert result.returncode == 0
    output = ".tact-pytest-" + stem + ".newick.tre"
    tacted = Tree.get(path=output, schema="newick", rooting="default-rooted")
    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)
    result = script_runner.run(
        "tact_check_results",
        output,
        "--taxonomy",
        taxonomy,
        "--backbone",
        backbone,
        "--output",
        ".tact-pytest-" + stem + ".check.csv",
        "--cores=1",
    )
    assert result.returncode == 0
    return (tacted, taxed, bbone)


@pytest.mark.parametrize("execution_number", execution_number)
def test_yule(script_runner, execution_number, datadir):
    backbone = os.path.join(datadir, "stem2.backbone.tre")
    taxonomy = os.path.join(datadir, "stem2.taxonomy.tre")
    taxed = Tree.get(path=taxonomy, schema="newick")
    bbone = Tree.get(path=backbone, schema="newick", rooting="default-rooted")
    result = script_runner.run(
        "tact_add_taxa",
        "--taxonomy",
        taxonomy,
        "--backbone",
        backbone,
        "--output",
        ".tact-pytest-yule",
        "-vv",
        "--yule",
    )
    assert result.returncode == 0
    output = ".tact-pytest-yule.newick.tre"
    tacted = Tree.get(path=output, schema="newick", rooting="default-rooted")
    ss = tacted.as_ascii_plot()
    sys.stderr.write(ss)
    result = script_runner.run(
        "tact_check_results",
        output,
        "--taxonomy",
        taxonomy,
        "--backbone",
        backbone,
        "--output",
        ".tact-pytest-yule.check.csv",
        "--cores=1",
    )
    assert result.returncode == 0
    return (tacted, taxed, bbone)


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "intrusion", "short_branch", "stem"])
def test_monophyly(script_runner, execution_number, datadir, stem):
    tacted, taxed, bbone = run_tact(script_runner, datadir, stem)
    extant = set([x.taxon.label for x in bbone.leaf_nodes()])
    for node in taxed.postorder_internal_node_iter(exclude_seed_node=True):
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


@pytest.mark.parametrize("execution_number", execution_number)
@pytest.mark.parametrize("stem", ["weirdness", "short_branch"])
def test_short_branch(script_runner, execution_number, datadir, stem):
    tacted, taxed, bbone = run_tact(script_runner, datadir, stem)
    n_short = 0
    for leaf in tacted.leaf_node_iter():
        if leaf.edge.length < 0.1:
            n_short += 1
    assert n_short <= 20


@pytest.mark.parametrize("execution_number", execution_number)
def test_stem_clade_attachment(script_runner, execution_number, datadir):
    tacted, taxed, bbone = run_tact(script_runner, datadir, "stem2")
    tacted.calc_node_ages()
    tacted.update_bipartitions()
    node = tacted.mrca(taxon_labels=["c1", "c2", "c3", "c4", "c5"])
    ages = [x.age < 15.16 for x in node.postorder_internal_node_iter()]
    assert any(ages)
