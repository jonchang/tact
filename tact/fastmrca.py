"""Fast MRCA lookups by taxon label."""

from collections.abc import Iterable

import dendropy

from .tree_util import get_tip_labels


class FastMRCA:
    """MRCA lookups by taxon label on one tree.

    `TaxonNamespace.taxa_bitmask(labels=...)` scans the whole namespace once per label, which
    is quadratic. This caches the bitmask of each label instead, and rebuilds the cache when taxa
    are added to or removed from the namespace.

    Args:
        tree: Tree to look up MRCAs in. Its bitmasks must be kept up to date, for example with
            `update_tree_view`, or by grafting with `graft_node(..., taxon_namespace=...)`.
    """

    def __init__(self, tree: dendropy.Tree) -> None:
        """See the class docstring for the arguments."""
        self.tree = tree
        # combined bitmask of the taxa carrying each label
        self._label_bitmasks: dict[str, int] = {}
        self._namespace_size = -1
        self._index_labels()

    def _index_labels(self) -> None:
        """Rebuild the label to bitmask index for the tree's taxon namespace."""
        tn = self.tree.taxon_namespace
        self._label_bitmasks.clear()
        for taxon in tn:
            if taxon.label is None:
                continue
            # normalize labels the same way DendroPy does when it matches them
            key = taxon.label if tn.is_case_sensitive else str(taxon.label).lower()
            self._label_bitmasks[key] = self._label_bitmasks.get(key, 0) | tn.taxon_bitmask(taxon)
        self._namespace_size = len(tn)

    def bitmask(self, labels: Iterable[str]) -> int:
        """Get a bitmask for the specified taxa labels.

        Intended to be equivalent to DendroPy's `TaxonNamespace.taxa_bitmask(labels=labels)`.

        * Labels that are not in the namespace are ignored
        * Case sensitivity of the namespace is respected

        Args:
            labels: Iterable of taxon labels.

        Returns:
            Bitmask representing the taxa.
        """
        tn = self.tree.taxon_namespace
        if len(tn) != self._namespace_size:
            # taxa were added to or removed from the namespace since the index was built
            self._index_labels()
        mask = 0
        if tn.is_case_sensitive:
            for label in labels:
                mask |= self._label_bitmasks.get(label, 0)
        else:
            for label in labels:
                mask |= self._label_bitmasks.get(str(label).lower(), 0)
        return mask

    def get(self, labels: Iterable[str]) -> dendropy.Node | None:
        """Get the MRCA node for the specified taxa.

        Returns the most recent common ancestor node if all descendants of that
        node are included in the label set, otherwise returns None.

        Args:
            labels: Iterable of taxon labels to find MRCA for.

        Returns:
            MRCA node if it forms a monophyletic group, None otherwise.
        """
        labels_set = set(labels)
        mrca = self.tree.mrca(leafset_bitmask=self.bitmask(labels_set))
        if mrca and labels_set.issuperset(get_tip_labels(mrca)):
            return mrca
        return None
