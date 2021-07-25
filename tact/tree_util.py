# -*- coding: utf-8 -*-

"""Functions specifically to handle DendroPy tree objects."""

import collections
import math

import dendropy

from .lib import optim_bd
from .lib import optim_yule


def get_birth_death_rates(node, sampfrac, yule=False, include_root=False):
    """
    Estimates the birth and death rates for the subtree descending from
    `node` with sampling fraction `sampfrac`. Optionally restrict to a
    Yule pure-birth model.
    """
    if yule:
        return optim_yule(get_ages(node, include_root), sampfrac)
    else:
        return optim_bd(get_ages(node, include_root), sampfrac)


def get_monophyletic_node(tree, species):
    """Returns the node or None that is the MRCA of the `species` in `tree`."""
    mrca = tree.mrca(taxon_labels=species)
    if not mrca:
        return None
    if mrca and species.issuperset(get_tip_labels(mrca)):
        return mrca


def get_ages(node, include_root=False):
    """Returns the list of ages of the children of a given `node`, optionally including the `node`'s age if `include_root` is True."""
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    if include_root:
        ages += [node.age]
    return ages


def get_tip_labels(tree_or_node):
    """Returns a `set` of tip labels for a node or tree."""
    try:
        return set([x.taxon.label for x in tree_or_node.leaf_node_iter()])
    except AttributeError:
        return set([x.taxon.label for x in tree_or_node.leaf_iter()])


def edge_iter(node, filter_fn=None):
    """
    Iterates over the child edge of `node` and all its descendants.
    Can optionally be filtered by `filter_fn`.
    """
    stack = list(node.child_edge_iter())
    while stack:
        edge = stack.pop()
        if filter_fn is None or filter_fn(edge):
            yield edge
        stack.extend(edge.head_node.child_edge_iter())


def get_tree(path, namespace=None):
    """
    Gets a DendroPy tree from a path and precalculate its node ages and bipartition bitmask.
    """
    tree = dendropy.Tree.get_from_path(
        path, schema="newick", taxon_namespace=namespace, rooting="default-rooted"
    )
    update_tree_view(tree)
    return tree


def update_tree_view(tree):
    """
    Mutates a DendroPy tree object with updated node ages and bipartition bitmask. We also
    correct for minor ultrametricity errors.

    Returns a list of tip labels.
    """
    tree.calc_node_ages(is_force_max_age=True)
    tree.update_bipartitions()
    return get_tip_labels(tree)


def is_binary(node):
    """Is the subtree under `node` a fully bifurcating tree?"""
    for x in node.preorder_internal_node_iter():
        if len(x.child_nodes()) != 2:
            return False
    return True


def is_ultrametric(tree, tolerance=1e-6):
    """Is the `tree` ultrametric, within a specified `tolerance`?

    Uses the relative difference between minimum and maximum root-to-tip distances.
    """
    tree.calc_node_root_distances()
    lengths = {}
    for leaf in tree.leaf_node_iter():
        lengths[leaf.taxon.label] = leaf.root_distance
    t_min = min(lengths.items(), key=lambda x: x[1])
    t_max = max(lengths.items(), key=lambda x: x[1])
    return (math.isclose(t_min[1], t_max[1], rel_tol=tolerance), (t_min, t_max))


def get_short_branches(node):
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge


def compute_node_depths(tree):
    res = dict()
    for leaf in tree.leaf_node_iter():
        cnt = 0
        for anc in leaf.ancestor_iter():
            if anc.label:
                cnt += 1
        res[leaf.taxon.label] = cnt
    return res


def ensure_tree_node_depths(tree):
    node_depths = compute_node_depths(tree)
    stats = collections.defaultdict(int)
    for v in node_depths.values():
        stats[v] += 1
    msg = ""
    if len(stats) > 1:
        msg += "The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:\n"
        for k in sorted(stats.keys()):
            msg += f"* {stats[k]} tips have {k} ranked ancestors\n"
    return msg
