"""Parity tests for the label index in `FastMRCA` with DendroPy."""

import os
import random

import dendropy
import pytest

from tact.fastmrca import FastMRCA
from tact.tree_util import get_monophyletic_node, get_tip_labels, update_tree_view

DATASETS = ["weirdness", "intrusion", "short_branch", "stem", "stem2"]


def dendropy_bitmask(tree, labels):
    """DendroPy reference implementation."""
    return tree.taxon_namespace.taxa_bitmask(labels=labels)


def load_pair(datadir, stem):
    """Load a taxonomy and backbone the way `tact_add_taxa` does, sharing one taxon namespace."""
    taxonomy = dendropy.Tree.get(
        path=os.path.join(datadir, stem + ".taxonomy.tre"), schema="newick", rooting="default-rooted"
    )
    backbone = dendropy.Tree.get(
        path=os.path.join(datadir, stem + ".backbone.tre"),
        schema="newick",
        rooting="default-rooted",
        taxon_namespace=taxonomy.taxon_namespace,
    )
    update_tree_view(backbone)
    return taxonomy, backbone


def label_sets(taxonomy, rng):
    """Every rank to be checked plus random subsets."""
    for node in taxonomy.preorder_internal_node_iter():
        yield get_tip_labels(node)
    all_labels = sorted(get_tip_labels(taxonomy))
    for _ in range(50):
        yield rng.sample(all_labels, rng.randint(1, len(all_labels)))


@pytest.mark.parametrize("stem", DATASETS)
def test_bitmask_matches_dendropy(datadir, stem):
    taxonomy, backbone = load_pair(datadir, stem)
    lookup = FastMRCA(backbone)
    for labels in label_sets(taxonomy, random.Random(stem)):
        assert lookup.bitmask(labels) == dendropy_bitmask(lookup.tree, labels)
        # any iterable of labels is accepted
        assert lookup.bitmask(iter(list(labels))) == dendropy_bitmask(lookup.tree, list(labels))


@pytest.mark.parametrize("stem", DATASETS)
def test_get_matches_reference(datadir, stem):
    taxonomy, backbone = load_pair(datadir, stem)
    lookup = FastMRCA(backbone)
    backbone_tips = get_tip_labels(backbone)
    checked = 0
    for labels in label_sets(taxonomy, random.Random(stem)):
        extant = backbone_tips.intersection(labels)
        if not extant:
            continue
        assert lookup.get(extant) is get_monophyletic_node(backbone, extant)
        checked += 1
    assert checked > 0


def test_bitmask_ignores_unknown_labels(datadir):
    _taxonomy, backbone = load_pair(datadir, "weirdness")
    lookup = FastMRCA(backbone)
    known = sorted(get_tip_labels(backbone))[:3]
    with_unknown = [*known, "Not a taxon"]
    assert lookup.bitmask(with_unknown) == dendropy_bitmask(lookup.tree, with_unknown) == lookup.bitmask(known)
    assert lookup.bitmask(["Not a taxon"]) == dendropy_bitmask(lookup.tree, ["Not a taxon"]) == 0
    assert lookup.bitmask([]) == 0


@pytest.mark.parametrize("is_case_sensitive", [False, True])
def test_bitmask_respects_namespace_case_sensitivity(is_case_sensitive):
    tn = dendropy.TaxonNamespace(is_case_sensitive=is_case_sensitive)
    tree = dendropy.Tree.get(
        data="((Aa_bb:1,Cc_dd:1):1,Ee_ff:2);",
        schema="newick",
        taxon_namespace=tn,
        case_sensitive_taxon_labels=is_case_sensitive,
    )
    update_tree_view(tree)
    lookup = FastMRCA(tree)
    for labels in (["Aa bb"], ["aa bb"], ["AA BB", "Ee ff"], ["cc DD", "Cc dd"]):
        assert lookup.bitmask(labels) == dendropy_bitmask(lookup.tree, labels)
    assert (lookup.bitmask(["aa bb"]) == 0) is is_case_sensitive


def test_bitmask_combines_taxa_sharing_a_label():
    tree = dendropy.Tree.get(data="((A:1,B:1):1,C:2);", schema="newick")
    tn = tree.taxon_namespace
    duplicate = tn.new_taxon("A")
    update_tree_view(tree)
    lookup = FastMRCA(tree)
    assert lookup.bitmask(["A"]) == dendropy_bitmask(lookup.tree, ["A"])
    assert lookup.bitmask(["A"]) & tn.taxon_bitmask(duplicate)


def test_bitmask_follows_namespace_changes():
    tree = dendropy.Tree.get(data="((A:1,B:1):1,C:2);", schema="newick")
    tn = tree.taxon_namespace
    update_tree_view(tree)
    lookup = FastMRCA(tree)
    assert lookup.bitmask(["D"]) == 0
    new_taxon = tn.new_taxon("D")
    assert lookup.bitmask(["D"]) == dendropy_bitmask(lookup.tree, ["D"]) == tn.taxon_bitmask(new_taxon)
    tn.remove_taxon(new_taxon)
    assert lookup.bitmask(["D"]) == dendropy_bitmask(lookup.tree, ["D"]) == 0


def assert_output_unchanged(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed):
    """`tact_add_taxa` must write the same files with the label index as with DendroPy's lookup."""
    indexed = add_taxa(backbone, taxonomy, tmp_path / "indexed", seed=seed)
    with monkeypatch.context() as patch:
        patch.setattr(FastMRCA, "bitmask", lambda self, labels: dendropy_bitmask(self.tree, labels))
        reference = add_taxa(backbone, taxonomy, tmp_path / "reference", seed=seed)
    assert indexed == reference
    # a different seed gives a different tree
    assert add_taxa(backbone, taxonomy, tmp_path / "other", seed=seed + 1)["newick.tre"] != indexed["newick.tre"]


@pytest.mark.parametrize("seed", range(3))
@pytest.mark.parametrize("stem", DATASETS)
def test_add_taxa_output_unchanged(monkeypatch, add_taxa, datadir, tmp_path, stem, seed):
    backbone = os.path.join(datadir, stem + ".backbone.tre")
    taxonomy = os.path.join(datadir, stem + ".taxonomy.tre")
    assert_output_unchanged(monkeypatch, add_taxa, tmp_path, backbone, taxonomy, seed)


def test_add_taxa_output_unchanged_carangaria(monkeypatch, add_taxa, carangaria, tmp_path):
    assert_output_unchanged(monkeypatch, add_taxa, tmp_path, *carangaria, 0)
