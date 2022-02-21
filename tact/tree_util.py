# -*- coding: utf-8 -*-

"""Functions specifically to handle DendroPy tree objects."""

import math
import random

import dendropy
import portion

from .lib import optim_bd
from .lib import optim_yule
from .exceptions import DisjointConstraintError


def get_birth_death_rates(node, sampfrac, yule=False, include_root=False):
    """
    Estimates the birth and death rates for the subtree descending from
    `node` with sampling fraction `sampfrac`. Optionally restrict to a
    Yule pure-birth model.
    """
    if yule:
        return optim_yule(get_ages(node, include_root), sampfrac)

    return optim_bd(get_ages(node, include_root), sampfrac)


def get_monophyletic_node(tree, species):
    """Returns the node or None that is the MRCA of the `species` in `tree`."""
    mrca = tree.mrca(taxon_labels=species)
    if mrca and species.issuperset(get_tip_labels(mrca)):
        return mrca

    return None


def get_ages(node, include_root=False):
    """
    Returns the list of ages of the children of a given `node`,
    optionally including the `node`'s age if `include_root` is True.
    """
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    if include_root:
        ages += [node.age]
    return ages


def get_tip_labels(tree_or_node):
    """Returns a `set` of tip labels for a node or tree."""
    try:
        return {x.taxon.label for x in tree_or_node.leaf_node_iter()}
    except AttributeError:
        return {x.taxon.label for x in tree_or_node.leaf_iter()}


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
    tree = dendropy.Tree.get_from_path(path, schema="newick", taxon_namespace=namespace, rooting="default-rooted")
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
    for internal_node in node.preorder_internal_node_iter():
        if len(internal_node.child_nodes()) != 2:
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
    """Yields an iterator of especially short edges under `node`."""
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge


def compute_node_depths(tree):
    """Returns a dictionary of node depths for each node with a label."""
    res = {}
    for leaf in tree.leaf_node_iter():
        cnt = 0
        for anc in leaf.ancestor_iter():
            if anc.label:
                cnt += 1
        res[leaf.taxon.label] = cnt
    return res


def graft_node(graft_recipient, graft, stem=False):
    """
    Grafts a node `graft` randomly in the subtree below node
    `graft_recipient`. The attribute `graft.age` must be set so
    we know where is the best place to graft the node. The node
    `graft` can optionally have child nodes, in this case the
    `edge.length` attribute should be set on all child nodes if
    the tree is to remain ultrametric.

    We graft things "below" a node by picking one of the children
    of that node and forcing it to be sister to the grafted node
    and adjusting the edge lengths accordingly. Therefore, the node
    *above* which the graft lives (i.e., the one that will be the child
    of the new graft) must fulfill the following requirements:

    1. Must not be the crown node (cannot graft things above crown node)
    2. Must be younger than the graft node (no negative branches)
    3. Seed node must be older than graft node (no negative branches)
    4. Must not be locked (intruding on monophyly)
    """

    def filter_fn(x):
        return x.head_node.age <= graft.age and x.head_node.parent_node.age >= graft.age and x.label != "locked"

    all_edges = list(edge_iter(graft_recipient))
    if stem:
        # also include the crown node's subtending edge
        all_edges.append(graft_recipient.edge)
    eligible_edges = [x for x in all_edges if filter_fn(x)]

    if not eligible_edges:
        raise Exception(f"could not place node {graft} in clade {graft_recipient}")

    focal_node = random.choice([x.head_node for x in eligible_edges])
    seed_node = focal_node.parent_node
    sisters = focal_node.sibling_nodes()

    # pick a child edge and detach its corresponding node
    #
    # DendroPy's Node.remove_child() messes with the edge lengths.
    # But, Node.clear_child_nodes() simply cuts that bit of the tree out.
    seed_node.clear_child_nodes()

    # set the correct edge length on the grafted node and make the grafted
    # node a child of the seed node
    graft.edge.length = seed_node.age - graft.age
    if graft.edge.length < 0:
        raise Exception("negative branch length")
    sisters.append(graft)
    seed_node.set_child_nodes(sisters)

    # make the focal node a child of the grafted node and set edge length
    focal_node.edge.length = graft.age - focal_node.age
    if focal_node.edge.length < 0:
        raise Exception("negative branch length")
    graft.add_child(focal_node)

    # return the (potentially new) crown of the clade
    if graft_recipient.parent_node == graft:
        return graft
    return graft_recipient


def lock_clade(node, stem=False):
    """
    Locks a clade descending from `node` so future grafts will avoid locked edges.
    """
    for edge in edge_iter(node):
        edge.label = "locked"
    if stem:
        node.edge.label = "locked"


def unlock_clade(node, stem=False):
    """
    Unlocks a clade descending from `node` so new tips can be grafted to its edges.
    """
    for edge in edge_iter(node):
        edge.label = ""
    if stem:
        node.edge.label = ""


def count_locked(node):
    """How many edges under `node` are locked?"""
    sum([x.label == "locked" for x in edge_iter(node)])


def is_fully_locked(node):
    """
    Are all the edges below `node` locked?
    """
    return all(x.label == "locked" for x in edge_iter(node))


def get_min_age(node):
    """
    Gets the minimum possible age that could be generated in a clade under `node`,
    assuming that grafts to locked edges are restricted.
    """
    interval = get_age_intervals(node)

    if not interval.atomic:
        raise DisjointConstraintError(f"Constraint on {node} implies disjoint interval {interval}")

    if interval.empty:
        return 0.0

    return interval.lower


def get_age_intervals(node):
    """
    Gets the (possibly disjoint) interval that could be generated in the
    clade under `node`, assuming that grafts to locked edges are restricted.
    """
    acc = portion.empty()
    for edge in edge_iter(node, lambda x: x.label != "locked"):
        acc = acc | portion.closed(edge.head_node.age, edge.tail_node.age)
    return acc
