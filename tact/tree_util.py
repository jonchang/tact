# -*- coding: utf-8 -*-

"""Functions specifically to handle DendroPy tree objects."""

import collections
import math
import random

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


def graft_node(graft_recipient, graft, stem=False):
    """
    Grafts a node `graft` randomly in the subtree below node
    `graft_recipient`. The attribute `graft.age` must be set so
    we know where is the best place to graft the node. The node
    `graft` can optionally have child nodes, in this case the
    `edge.length` attribute should be set on all child nodes if
    the tree is to remain ultrametric.
    """

    # We graft things "below" a node by picking one of the children
    # of that node and forcing it to be sister to the grafted node
    # and adjusting the edge lengths accordingly. Therefore, the node
    # *above* which the graft lives (i.e., the one that will be the child
    # of the new graft) must fulfill the following requirements:
    #
    # 1. Must not be the crown node (cannot graft things above crown node)
    # 2. Must be younger than the graft node (no negative branches)
    # 3. Seed node must be older than graft node (no negative branches)
    # 4. Must not be locked (intruding on monophyly)
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


def lock_clade(node):
    """
    "Locks a clade descending from `node` so future grafts will avoid locked edges.
    """
    for edge in edge_iter(node):
        edge.label = "locked"


def unlock_clade(node):
    for edge in edge_iter(node):
        edge.label = ""


def count_locked(node):
    sum([x.label == "locked" for x in edge_iter(node)])


def is_fully_locked(node):
    return all([x.label == "locked" for x in edge_iter(node)])
