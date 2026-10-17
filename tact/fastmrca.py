"""Singleton object that helps speed up MRCA lookups."""

from typing import TYPE_CHECKING

from .tree_util import get_tip_labels

if TYPE_CHECKING:
    import dendropy

tree: "dendropy.Tree | None" = None

# Maps each label in a DendroPy tree to the combined bitmask of the taxa carrying that label.
# `TaxonNamespace.taxa_bitmask(labels=...)` scans the whole namespace once per label, which
# is quadratic. Cache per-label bitmasks to improve lookup speed for large clades.
_label_bitmasks: dict[str, int] = {}
_namespace_size = -1


def _index_labels(tn: "dendropy.TaxonNamespace") -> None:
    """Rebuild the label to bitmask index for a taxon namespace."""
    global _namespace_size
    _label_bitmasks.clear()
    for taxon in tn:
        if taxon.label is None:
            continue
        # normalize labels the same way DendroPy does when it matches them
        key = taxon.label if tn.is_case_sensitive else str(taxon.label).lower()
        _label_bitmasks[key] = _label_bitmasks.get(key, 0) | tn.taxon_bitmask(taxon)
    _namespace_size = len(tn)


def initialize(phy: "dendropy.Tree") -> None:
    """Initialize the fastmrca singleton with a tree.

    Args:
        phy: DendroPy tree object to use for MRCA lookups.
    """
    global tree
    tree = phy
    _index_labels(phy.taxon_namespace)


def bitmask(labels: list[str]) -> int:
    """Get a bitmask for the specified taxa labels.

    Intended to be equivalent to DendroPy's `TaxonNamespace.taxa_bitmask(labels=labels)`.

    * Labels that are not in the namespace are ignored
    * Case sensitivity of the namespace is respected

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
    if len(tn) != _namespace_size:
        # taxa were added to or removed from the namespace since the index was built
        _index_labels(tn)
    mask = 0
    if tn.is_case_sensitive:
        for label in labels:
            mask |= _label_bitmasks.get(label, 0)
    else:
        for label in labels:
            mask |= _label_bitmasks.get(str(label).lower(), 0)
    return mask


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
