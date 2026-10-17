"""Checks that grafting keeps the tree ready for MRCA lookups without refreshing the whole tree."""

import os
import random

import dendropy
import pytest

from tact import cli_add_taxa, fastmrca
from tact.tree_util import get_tip_labels, graft_node, update_tree_view

DATASETS = ["weirdness", "intrusion", "short_branch", "stem", "stem2"]


def assert_tree_view_current(tree):
    """Assert that stored bitmasks and ages are what `update_tree_view` would compute."""
    tn = tree.taxon_namespace
    expected = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            expected[node] = tn.taxon_bitmask(node.taxon)
            assert node.age == 0
        else:
            expected[node] = 0
            for child in node.child_node_iter():
                expected[node] |= expected[child]
            # backbones are only ultrametric to a tolerance, so ages are taken from the longest path
            oldest = max(child.age + child.edge.length for child in node.child_node_iter())
            assert node.age == pytest.approx(oldest, abs=1e-9)
        assert node.edge.bipartition.leafset_bitmask == expected[node]


def new_cherry(tn, label, age):
    """A new internal node of a given age subtending one new tip, as built by `fill_new_taxa`."""
    node = dendropy.Node()
    node.age = age
    leaf = node.new_child(taxon=tn.require_taxon(label), edge_length=age)
    leaf.age = 0
    return node


@pytest.fixture
def small_tree():
    tn = dendropy.TaxonNamespace(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"])
    tree = dendropy.Tree.get(
        data="((A:1,B:1):2,(C:2,D:2):1);", schema="newick", rooting="default-rooted", taxon_namespace=tn
    )
    update_tree_view(tree)
    return tree


def mrca_labels(tree, labels):
    return get_tip_labels(tree.mrca(leafset_bitmask=tree.taxon_namespace.taxa_bitmask(labels=labels)))


def test_graft_tips_keeps_tree_view_current(small_tree):
    tn = small_tree.taxon_namespace
    random.seed(0)
    root = small_tree.seed_node
    for label, age in [("E", 0.5), ("F", 1.5), ("G", 2.5), ("H", 0.1), ("I", 2.9)]:
        graft_node(root, new_cherry(tn, label, age), taxon_namespace=tn)
        assert_tree_view_current(small_tree)
        assert label in get_tip_labels(small_tree)
    assert mrca_labels(small_tree, list("ABCDEFGHI")) == set("ABCDEFGHI")


def test_graft_on_stem_keeps_tree_view_current(small_tree):
    tn = small_tree.taxon_namespace
    cherry = small_tree.mrca(taxon_labels=["A", "B"])
    # older than the crown of (A, B), so the only eligible edge is its stem
    crown = graft_node(cherry, new_cherry(tn, "E", 1.5), stem=True, taxon_namespace=tn)
    assert crown is not cherry
    assert_tree_view_current(small_tree)
    assert mrca_labels(small_tree, ["A", "E"]) == {"A", "B", "E"}
    assert mrca_labels(small_tree, ["A", "B"]) == {"A", "B"}


def test_graft_clade_keeps_tree_view_current(small_tree):
    tn = small_tree.taxon_namespace
    random.seed(0)
    clade = cli_add_taxa.create_clade(tn, ["E", "F", "G", "H"], [2.5, 2.0, 1.0, 0.5])
    graft_node(small_tree.seed_node, clade.seed_node, taxon_namespace=tn)
    assert_tree_view_current(small_tree)
    assert mrca_labels(small_tree, ["E", "F", "G", "H"]) == {"E", "F", "G", "H"}
    # a second graft on a tree whose bitmasks were only ever updated incrementally
    graft_node(small_tree.seed_node, new_cherry(tn, "I", 0.25), taxon_namespace=tn)
    assert_tree_view_current(small_tree)


def test_graft_without_namespace_needs_full_refresh(small_tree):
    tn = small_tree.taxon_namespace
    random.seed(0)
    before = small_tree.seed_node.edge.bipartition.leafset_bitmask
    graft_node(small_tree.seed_node, new_cherry(tn, "E", 0.5))
    assert small_tree.seed_node.edge.bipartition.leafset_bitmask == before
    update_tree_view(small_tree)
    assert_tree_view_current(small_tree)


def clade_ages(newick, tn):
    """Map every clade in a tree, as a leafset bitmask, to its age."""
    tree = dendropy.Tree.get(data=newick.decode(), schema="newick", rooting="default-rooted", taxon_namespace=tn)
    update_tree_view(tree)
    return {node.edge.bipartition.leafset_bitmask: node.age for node in tree.preorder_node_iter()}


def check_every_graft(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed):
    """The tree view must be current after each graft over a whole run of `tact_add_taxa`."""
    grafts = 0

    def checked_graft_node(*args, **kwargs):
        nonlocal grafts
        result = graft_node(*args, **kwargs)
        assert_tree_view_current(fastmrca.tree)
        grafts += 1
        return result

    monkeypatch.setattr(cli_add_taxa, "graft_node", checked_graft_node)
    add_taxa(backbone, taxonomy, tmp_path / "checked", seed=seed)
    assert grafts > 0


def check_matches_full_refresh(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed):
    """`tact_add_taxa` must build the same tree as it did when refreshing the whole tree per graft."""

    def refreshing_graft_node(*args, **kwargs):
        result = graft_node(*args, **kwargs)
        update_tree_view(fastmrca.tree)
        return result

    incremental = add_taxa(backbone, taxonomy, tmp_path / "incremental", seed=seed)
    with monkeypatch.context() as patch:
        patch.setattr(cli_add_taxa, "graft_node", refreshing_graft_node)
        reference = add_taxa(backbone, taxonomy, tmp_path / "reference", seed=seed)

    assert incremental["rates.csv"] == reference["rates.csv"]
    tn = dendropy.TaxonNamespace()
    ours = clade_ages(incremental["newick.tre"], tn)
    theirs = clade_ages(reference["newick.tre"], tn)
    # same clades, therefore same topology and tips
    assert ours.keys() == theirs.keys()
    # recomputing ages from branch lengths after each graft only moves them by rounding error
    for clade, age in ours.items():
        assert age == pytest.approx(theirs[clade], abs=1e-9)
    # a different seed gives a different tree
    assert add_taxa(backbone, taxonomy, tmp_path / "other", seed=seed + 1)["newick.tre"] != incremental["newick.tre"]


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("stem", DATASETS)
def test_add_taxa_tree_view_current_after_every_graft(monkeypatch, add_taxa, datadir, tmp_path, stem, seed):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    check_every_graft(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed)


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("stem", DATASETS)
def test_add_taxa_matches_full_refresh(monkeypatch, add_taxa, datadir, tmp_path, stem, seed):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    check_matches_full_refresh(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed)


def test_add_taxa_tree_view_current_after_every_graft_carangaria(monkeypatch, add_taxa, carangaria, tmp_path):
    check_every_graft(monkeypatch, add_taxa, tmp_path, *carangaria, 0)


def test_add_taxa_matches_full_refresh_carangaria(monkeypatch, add_taxa, carangaria, tmp_path):
    check_matches_full_refresh(monkeypatch, add_taxa, tmp_path, *carangaria, 0)
