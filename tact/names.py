"""Tooling to reconcile species and taxon names.

TACT shows names with spaces and writes them to files with underscores.
Two names that differ only in spaces, underscores or case are the same name.
`name_key` gives the key they share.
"""

import re
from collections.abc import Iterable, Mapping

import dendropy

# Runs of whitespace (including the non-breaking spaces spreadsheets add) and underscores
_SEPARATORS = re.compile(r"[\s_]+")


def clean_name(name: str) -> str:
    """Spell `name` with single spaces between words.

    Intended to be used for output written to the screen.
    """
    return _SEPARATORS.sub(" ", name).strip()


def name_key(name: str) -> str:
    """Key shared by every spelling of `name` that differs only in spaces, underscores or case."""
    return clean_name(name).casefold()


def written_name(name: str) -> str:
    """Spell a clean `name` with underscores.

    Intended to be used for output written to files.
    """
    return name.replace(" ", "_")


def spelling_clashes(names: Iterable[str]) -> list[list[str]]:
    """Group the distinct `names` that share a `name_key`, keeping the groups of more than one."""
    groups: dict[str, list[str]] = {}
    for name in names:
        group = groups.setdefault(name_key(name), [])
        if name not in group:
            group.append(name)
    return [group for group in groups.values() if len(group) > 1]


def clean_tree_names(tree: dendropy.Tree, keep: Iterable[str] = ()) -> tuple[list[list[str]], list[list[str]]]:
    """Clean the tip and internal node labels of `tree` in place, except the internal labels in `keep`.

    Returns:
        Groups of tip labels, and of internal node labels, that were spelled differently but share a `name_key`.
    """
    keep = set(keep)
    tips = [leaf.taxon for leaf in tree.leaf_node_iter()]
    nodes = [node for node in tree.preorder_internal_node_iter() if node.label and node.label not in keep]
    clashes = (spelling_clashes(taxon.label for taxon in tips), spelling_clashes(node.label for node in nodes))
    for taxon in tips:
        taxon.label = clean_name(taxon.label)
    for node in nodes:
        node.label = clean_name(node.label)
    return clashes


def adopt_names(
    tree: dendropy.Tree, namespace: dendropy.TaxonNamespace, add_missing: bool = False
) -> tuple[list[str], dict[str, str]]:
    """Move the tips of `tree` onto the taxa of `namespace` that have the same `name_key`.

    Args:
        tree: A tree whose labels have been cleaned.
        namespace: The taxa to use.
        add_missing: Add a taxon to `namespace` for each tip that has none, instead of leaving the tip as it was.

    Returns:
        The labels of the tips with no taxon in `namespace`, and the labels spelled differently in `namespace`,
        mapped to their spelling there.
    """
    taxa: Mapping[str, dendropy.Taxon] = {name_key(taxon.label): taxon for taxon in namespace}
    missing = []
    respelled = {}
    for leaf in tree.leaf_node_iter():
        label = leaf.taxon.label
        taxon = taxa.get(name_key(label))
        if taxon is None:
            missing.append(label)
            if not add_missing:
                continue
            taxon = namespace.new_taxon(label)
        elif taxon.label != label:
            respelled[label] = taxon.label
        leaf.taxon = taxon
    if add_missing or not missing:
        tree.taxon_namespace = namespace
    return missing, respelled
