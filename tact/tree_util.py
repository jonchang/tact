"""Functions specifically to handle DendroPy tree objects."""

import math
import random
import re
from collections.abc import Callable, Iterator

import dendropy
import portion

from .exceptions import DisjointConstraintError
from .lib import optim_bd, optim_yule


def get_birth_death_rates(
    node: dendropy.Node, sampfrac: float, yule: bool = False, include_root: bool = False
) -> tuple[float, float]:
    """Estimate birth-death rates from a subtree.

    Computes maximum likelihood birth and death rates for the subtree
    descending from the given node, optionally using a Yule (pure birth) model.

    Args:
        node: DendroPy node representing the root of the subtree.
        sampfrac: Sampling fraction in (0, 1].
        yule: If True, use Yule model (zero extinction). Default: False.
        include_root: If True, include root node age in calculations. Default: False.

    Returns:
        Tuple of (birth rate, death rate).
    """
    if yule:
        return optim_yule(get_ages(node, include_root), sampfrac)

    return optim_bd(get_ages(node, include_root), sampfrac)


def get_monophyletic_node(tree: dendropy.Tree, species: set[str]) -> dendropy.Node | None:
    """Get the MRCA node if it forms a monophyletic group.

    Args:
        tree: DendroPy tree object.
        species: Set of taxon labels to check for monophyly.

    Returns:
        MRCA node if all its descendants are in the species set, None otherwise.
    """
    mrca = tree.mrca(taxon_labels=species)
    if mrca and species.issuperset(get_tip_labels(mrca)):
        return mrca

    return None


def get_ages(node: dendropy.Node, include_root: bool = False) -> list[float]:
    """Get list of node ages in a subtree.

    Args:
        node: DendroPy node representing the root of the subtree.
        include_root: If True, include the root node's age. Default: False.

    Returns:
        List of node ages, sorted in descending order.
    """
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    if include_root:
        ages += [node.age]
    return ages


def get_tip_labels(tree_or_node: dendropy.Tree | dendropy.Node) -> set[str]:
    """Get set of tip labels for a tree or node.

    Args:
        tree_or_node: DendroPy tree or node object.

    Returns:
        Set of taxon labels for all tips.
    """
    if isinstance(tree_or_node, dendropy.Tree):
        return {x.taxon.label for x in tree_or_node.leaf_node_iter()}
    return {x.taxon.label for x in tree_or_node.leaf_iter()}


def edge_iter(
    node: dendropy.Node, filter_fn: Callable[[dendropy.Edge], bool] | None = None
) -> Iterator[dendropy.Edge]:
    """Iterate over edges in a subtree.

    Args:
        node: DendroPy node representing the root of the subtree.
        filter_fn: Optional function to filter edges. Default: None.

    Yields:
        DendroPy edge objects in the subtree.
    """
    stack = list(node.child_edge_iter())
    while stack:
        edge = stack.pop()
        if filter_fn is None or filter_fn(edge):
            yield edge
        stack.extend(edge.head_node.child_edge_iter())


# Quoted labels (with '' for a quote) and [comments], which may contain semicolons
_NEWICK_QUOTES_AND_COMMENTS = re.compile(r"'(?:[^']|'')*'|\[[^\]]*\]")


def count_newick_trees(path: str) -> int:
    """Count the trees in a Newick file by their terminating semicolons.

    Args:
        path: File path to the tree file.

    Returns:
        Number of trees in the file.
    """
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    return _NEWICK_QUOTES_AND_COMMENTS.sub("", text).count(";")


def get_tree(path: str, namespace: dendropy.TaxonNamespace | None = None) -> dendropy.Tree:
    """Load a DendroPy tree from a file path.

    Loads a Newick-format tree and precalculates node ages and bipartition
    bitmasks for efficient operations.

    Args:
        path: File path to the tree file.
        namespace: Optional taxon namespace to use. Default: None.

    Returns:
        DendroPy tree object with ages and bipartitions calculated.
    """
    tree = dendropy.Tree.get_from_path(path, schema="newick", taxon_namespace=namespace, rooting="default-rooted")
    update_tree_view(tree)
    return tree


def update_tree_view(tree: dendropy.Tree) -> set[str]:
    """Update a tree with node ages and bipartition bitmasks.

    Performs in-place updates to calculate node ages and bipartition bitmasks,
    correcting for minor ultrametricity errors.

    Args:
        tree: DendroPy tree object to update.

    Returns:
        Set of tip labels in the tree.
    """
    tree.calc_node_ages(is_force_max_age=True)
    tree.update_bipartitions()
    return get_tip_labels(tree)


def is_binary(node: dendropy.Node) -> bool:
    """Check if a subtree is fully bifurcating.

    Args:
        node: DendroPy node representing the root of the subtree.

    Returns:
        True if all internal nodes have exactly two children, False otherwise.
    """
    for internal_node in node.preorder_internal_node_iter():
        if len(internal_node.child_nodes()) != 2:
            return False
    return True


def is_ultrametric(
    tree: dendropy.Tree, tolerance: float = 1e-6
) -> tuple[bool, tuple[tuple[str, float], tuple[str, float]]]:
    """Check if a tree is ultrametric within a specified tolerance.

    Args:
        tree: DendroPy tree object to check.
        tolerance: Relative tolerance for ultrametricity check. Default: 1e-6.

    Returns:
        Tuple of (is_ultrametric, (min_tip, max_tip)) where min_tip and max_tip
        are tuples of (label, root_distance) for the tips with minimum and
        maximum root-to-tip distances.
    """
    tree.calc_node_root_distances()
    lengths: dict[str, float] = {}
    for leaf in tree.leaf_node_iter():
        lengths[leaf.taxon.label] = leaf.root_distance
    t_min = min(lengths.items(), key=lambda x: x[1])
    t_max = max(lengths.items(), key=lambda x: x[1])
    return (math.isclose(t_min[1], t_max[1], rel_tol=tolerance), (t_min, t_max))


def get_short_branches(node: dendropy.Node) -> Iterator[dendropy.Edge]:
    """Get edges with very short branch lengths.

    Args:
        node: DendroPy node representing the root of the subtree.

    Yields:
        DendroPy edge objects with length <= 0.001.
    """
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge


def compute_node_depths(tree: dendropy.Tree) -> dict[str, int]:
    """Compute node depths for all labeled nodes.

    Depth is defined as the number of labeled ancestor nodes.

    Args:
        tree: DendroPy tree object.

    Returns:
        Dictionary mapping tip labels to their node depths.
    """
    res: dict[str, int] = {}
    for leaf in tree.leaf_node_iter():
        cnt = 0
        for anc in leaf.ancestor_iter():
            if anc.label:
                cnt += 1
        res[leaf.taxon.label] = cnt
    return res


def set_leafset_bitmask(node: dendropy.Node, bitmask: int) -> None:
    """Set the leafset bitmask of the edge subtending a node.

    Args:
        node: Node whose subtending edge should be updated.
        bitmask: Bitmask of the taxa descending from the node.
    """
    node.edge.bipartition = dendropy.Bipartition(
        leafset_bitmask=bitmask, compile_bipartition=False, is_mutable=True, is_rooted=True
    )


def encode_leafset_bitmasks(node: dendropy.Node, taxon_namespace: dendropy.TaxonNamespace) -> int:
    """Compute and store leafset bitmasks for every edge in a subtree.

    Args:
        node: Node representing the root of the subtree.
        taxon_namespace: Taxon namespace that the leaves of the subtree belong to.

    Returns:
        Bitmask of the taxa in the subtree.
    """
    bitmask = 0
    for nd in node.postorder_iter():
        if nd.is_leaf():
            bitmask = taxon_namespace.taxon_bitmask(nd.taxon)
        else:
            bitmask = 0
            for child in nd.child_node_iter():
                bitmask |= child.edge.bipartition.leafset_bitmask
        set_leafset_bitmask(nd, bitmask)
    return bitmask


def graft_node(
    graft_recipient: dendropy.Node,
    graft: dendropy.Node,
    stem: bool = False,
    taxon_namespace: dendropy.TaxonNamespace | None = None,
) -> dendropy.Node:
    """Graft a node randomly into a subtree.

    Grafts a node into the subtree below the recipient node, randomly selecting
    an eligible edge. The graft node's age must be set. Edge lengths are
    adjusted to maintain ultrametricity.

    The eligible edge (above which the graft is placed) must satisfy:

    1. Not be the crown node
    2. Head node younger than graft node
    3. Tail node older than graft node
    4. Not be locked

    If `taxon_namespace` is given, ensure that the tree is left ready for MRCA lookups
    without a subsequent call to `update_tree_view`. Prior to grafting, the bitmasks of
    the tree must be up-to-date, and every node in the grafted subtree must already have
    its age set.

    * leafset bitmasks are computed for the grafted subtree
    * new taxa are added to the bitmask of every edge between the graft and the root

    Args:
        graft_recipient: Node representing the clade to graft into.
        graft: Node to graft (must have age attribute set).
        stem: If True, allow grafting on the stem edge. Default: False.
        taxon_namespace: If given, keep the leafset bitmasks of the tree up to date. Default: None.

    Returns:
        The crown node of the clade (may be the graft if it becomes the new crown).

    Raises:
        Exception: If no eligible edge is found or negative branch lengths result.
    """

    def filter_fn(x: dendropy.Edge) -> bool:
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

    if taxon_namespace is not None:
        new_taxa = 0
        for child in graft.child_node_iter():
            if child is not focal_node:
                new_taxa |= encode_leafset_bitmasks(child, taxon_namespace)
        set_leafset_bitmask(graft, new_taxa | focal_node.edge.bipartition.leafset_bitmask)
        for ancestor in graft.ancestor_iter():
            set_leafset_bitmask(ancestor, new_taxa | ancestor.edge.bipartition.leafset_bitmask)

    # return the (potentially new) crown of the clade
    if graft_recipient.parent_node == graft:
        return graft
    return graft_recipient


def lock_clade(node: dendropy.Node, stem: bool = False) -> None:
    """Lock a clade to prevent future grafts.

    Marks all edges in the subtree as locked, preventing new grafts from
    being placed on them.

    Args:
        node: Node representing the root of the clade to lock.
        stem: If True, also lock the stem edge. Default: False.
    """
    for edge in edge_iter(node):
        edge.label = "locked"
    if stem:
        node.edge.label = "locked"


def unlock_clade(node: dendropy.Node, stem: bool = False) -> None:
    """Unlock a clade to allow future grafts.

    Removes the locked label from all edges in the subtree.

    Args:
        node: Node representing the root of the clade to unlock.
        stem: If True, also unlock the stem edge. Default: False.
    """
    for edge in edge_iter(node):
        edge.label = ""
    if stem:
        node.edge.label = ""


def is_fully_locked(node: dendropy.Node) -> bool:
    """Check if all edges in a subtree are locked.

    Args:
        node: Node representing the root of the subtree.

    Returns:
        True if all edges are locked, False otherwise.
    """
    return all(x.label == "locked" for x in edge_iter(node))


def node_age(node: dendropy.Node) -> float:
    """Get the age of a node.

    DendroPy initializes `Node.age` to None, so this narrows it to a float for type checkers.

    Args:
        node: DendroPy node with its age set.

    Returns:
        The age of the node.

    Raises:
        ValueError: If the node's age has not been set.
    """
    if node.age is None:
        raise ValueError(f"{node} has no age; have node ages been calculated?")
    return node.age


def get_min_age(node: dendropy.Node) -> float:
    """Get the minimum possible age for a graft in a clade.

    Computes the minimum age that could be generated by grafting into the
    clade, considering only unlocked edges.

    Args:
        node: Node representing the root of the clade.

    Returns:
        Minimum possible age (0.0 if interval is empty).

    Raises:
        DisjointConstraintError: If the age interval is not atomic (disjoint).
    """
    interval = get_age_intervals(node)

    if not interval.atomic:
        raise DisjointConstraintError(f"Constraint on {node} implies disjoint interval {interval}")

    if interval.empty:
        return 0.0

    return float(interval.lower)


def get_age_intervals(node: dendropy.Node) -> portion.Interval:
    """Get the age interval for possible grafts in a clade.

    Computes the union of all age intervals from unlocked edges, which
    represents all possible ages where a graft could be placed.

    Args:
        node: Node representing the root of the clade.

    Returns:
        Portion interval (possibly disjoint) representing valid graft ages.
    """
    acc = portion.empty()
    for edge in edge_iter(node, lambda x: x.label != "locked"):
        acc = acc | portion.closed(edge.head_node.age, edge.tail_node.age)
    return acc
