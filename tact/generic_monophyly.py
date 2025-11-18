"""Module that provides a generic monophyly singleton data object."""

import re
from collections import defaultdict
from typing import TYPE_CHECKING

from .tree_util import get_monophyletic_node, get_tip_labels

if TYPE_CHECKING:
    import dendropy

_valid_monophyly: dict[str, set[str]] = {}


def init(tree: "dendropy.Tree") -> dict[str, set[str]]:
    """Initialize a monophyly singleton object.

    Builds a dictionary mapping genus names to sets of species that form
    monophyletic groups in the tree. Genus names are extracted from tip
    labels by splitting on underscores or spaces.

    Args:
        tree: DendroPy tree object to analyze.

    Returns:
        Dictionary mapping genus names to sets of species labels that
        form monophyletic groups.
    """
    if len(_valid_monophyly):
        return _valid_monophyly

    _genera_map: defaultdict[str, list[str]] = defaultdict(list)
    tips = get_tip_labels(tree)
    for tip in tips:
        genus, _ = re.split("[_ ]+", tip, maxsplit=1)
        _genera_map[genus].append(tip)

    for genus, species in _genera_map.items():
        species_set = set(species)
        node = get_monophyletic_node(tree, species_set)
        if node:
            _valid_monophyly[genus] = species_set
    return _valid_monophyly
