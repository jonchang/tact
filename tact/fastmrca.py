"""Singleton object that helps speed up MRCA lookups."""

from typing import TYPE_CHECKING

from .tree_util import get_tip_labels

if TYPE_CHECKING:
    import dendropy

tree: "dendropy.Tree | None" = None


def initialize(phy: "dendropy.Tree") -> None:
    """Initialize the fastmrca singleton with a tree.

    Args:
        phy: DendroPy tree object to use for MRCA lookups.
    """
    global tree
    tree = phy


def bitmask(labels: list[str]) -> int:
    """Get a bitmask for the specified taxa labels.

    Args:
        labels: List of taxon labels.

    Returns:
        Bitmask representing the taxa.

    Raises:
        RuntimeError: If fastmrca has not been initialized.
    """
    global tree
    if tree is None:
        raise RuntimeError("fastmrca not initialized")
    tn = tree.taxon_namespace
    return tn.taxa_bitmask(labels=labels)  # type: ignore[no-any-return]


def get(labels: list[str]) -> "dendropy.Node | None":
    """Get the MRCA node for the specified taxa.

    Returns the most recent common ancestor node if all descendants of that
    node are included in the label set, otherwise returns None.

    Args:
        labels: List of taxon labels to find MRCA for.

    Returns:
        MRCA node if it forms a monophyletic group, None otherwise.

    Raises:
        RuntimeError: If fastmrca has not been initialized.
    """
    global tree
    if tree is None:
        raise RuntimeError("fastmrca not initialized")
    labels_set = set(labels)
    mrca = tree.mrca(leafset_bitmask=bitmask(labels))
    if mrca and labels_set.issuperset(get_tip_labels(mrca)):
        return mrca
    return None


def fastmrca_getter(tn: "dendropy.TaxonNamespace", x: list[str]) -> int:
    """Helper function to compute bitmask for parallel processing.

    Args:
        tn: Taxon namespace object.
        x: List of taxon labels.

    Returns:
        Bitmask representing the taxa.
    """
    taxa = tn.get_taxa(labels=x)
    mask = 0
    for taxon in taxa:
        mask |= int(tn.taxon_bitmask(taxon))
    return mask
