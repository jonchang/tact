import pytest
import os

from tact.tree_util import get_age_intervals, get_tree, lock_clade, unlock_clade

def test_disjoint(datadir):
    disjoint_tree = os.path.join(datadir, "disjoint.tre")
    tree = get_tree(disjoint_tree)

    lock_clade(tree.seed_node)
    unlock_clade(tree.mrca(taxon_labels=["A1", "A2"]))
    unlock_clade(tree.mrca(taxon_labels=["B1", "B2", "C1", "C2"]))
    lock_clade(tree.mrca(taxon_labels=["B1", "B2"]))
    lock_clade(tree.mrca(taxon_labels=["C1", "C2"]))

    res = get_age_intervals(tree.seed_node)
    assert not res.atomic
