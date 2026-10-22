"""Record where tact places each species and which taxa it had to skip."""

import csv
from collections.abc import Iterable
from typing import Self

import dendropy

FIELDS = (
    "species",
    "expected_taxon",
    "placed_taxon",
    "ranks_up",
    "new_clade",
    "skipped_taxa",
    "n_intruders",
    "intruders",
)

# max intruders listed per row of the CSV
MAX_INTRUDERS = 5

# taxa listed in the summary
MAX_SUMMARY_TAXA = 10


def fmt_species_list(spp: Iterable[str]) -> str:
    """Name the first two species and count the rest."""
    spp = list(spp)
    if len(spp) == 3:
        return f"{spp[0]}, {spp[1]} and 1 other"
    if len(spp) > 3:
        return f"{spp[0]}, {spp[1]} and {len(spp) - 2} others"
    return " and ".join(spp)


class PlacementReport:
    """Writes a CSV row per placed species, and summarizes placements at the end.

    A species is expected to be placed in the lowest labeled taxon containing it
    that has species in the backbone. When that taxon is not monophyletic in the
    backbone, missing species are instead added to the first ancestor that is.
    Record the skipped taxa and the backbone tips that break up the expected taxon.

    Args:
        path: CSV file to write.
        taxonomy: Taxonomy tree.
        backbone_tips: Labels of the tips in the backbone before any species are added.
    """

    def __init__(self, path: str, taxonomy: dendropy.Tree, backbone_tips: Iterable[str]) -> None:
        """See the class docstring for the arguments."""
        self.path = path
        self.backbone_tips = frozenset(backbone_tips)
        self.added = 0
        self.not_monophyletic: dict[str, list[str]] = {}
        # expected taxon => [number of species, taxon they were placed in, ranks up]
        self.fallbacks: dict[str, list] = {}
        self._parents: dict[str, str | None] = {}
        self._expected: dict[str, str] = {}
        self._index_taxonomy(taxonomy)
        self._file = open(path, "w", newline="", encoding="utf-8")
        self._writer = csv.writer(self._file)
        self._writer.writerow(FIELDS)
        self._file.flush()

    def _index_taxonomy(self, taxonomy: dendropy.Tree) -> None:
        """Find the expected taxon of each species missing from the backbone, and each taxon's labeled parent."""
        sampled = {}
        for node in taxonomy.postorder_node_iter():
            if node.is_leaf():
                sampled[node] = node.taxon.label in self.backbone_tips
            else:
                sampled[node] = any(sampled[child] for child in node.child_node_iter())
        for node in taxonomy.preorder_internal_node_iter():
            if node.label:
                parent = node.parent_node
                while parent is not None and not parent.label:
                    parent = parent.parent_node
                self._parents[node.label] = parent.label if parent is not None else None
        for leaf in taxonomy.leaf_node_iter():
            if sampled[leaf]:
                continue
            node = leaf.parent_node
            while node is not None and not (node.label and sampled[node]):
                node = node.parent_node
            if node is not None:
                self._expected[leaf.taxon.label] = node.label

    def __enter__(self) -> Self:
        """Return the report."""
        return self

    def __exit__(self, *exc_info: object) -> None:
        """Close the CSV file."""
        self.close()

    def close(self) -> None:
        """Close the CSV file."""
        self._file.close()

    def taxon_not_monophyletic(self, taxon: str, mrca_tips: Iterable[str], species: Iterable[str]) -> list[str]:
        """Record a skipped taxon because it is not monophyletic.

        Args:
            taxon: Label of the taxon.
            mrca_tips: Tips below the MRCA of the taxon's species in the tree.
            species: All species in the taxon.

        Returns:
            Sorted labels of the backbone tips below the MRCA that do not belong to the taxon.
        """
        intruders = sorted((set(mrca_tips) & self.backbone_tips) - set(species))
        self.not_monophyletic[taxon] = intruders
        return intruders

    def placed(self, species: Iterable[str], taxon: str, new_clade: str | None = None) -> None:
        """Record and write out species just added to the tree.

        Args:
            species: Labels of the added species.
            taxon: Taxon whose clade in the tree they were added to.
            new_clade: Label of the clade generated for them, if they were added as a new clade.

        Raises:
            ValueError: If `taxon` does not contain the expected taxon of a species.
        """
        for name in sorted(species):
            expected = self._expected.get(name, taxon)
            skipped = []
            ranks_up = 0
            current = expected
            while current != taxon:
                if current is None:
                    raise ValueError(
                        f"{name} was placed in {taxon}, which does not contain its expected taxon {expected}"
                    )
                if current in self.not_monophyletic:
                    skipped.append(current)
                current = self._parents.get(current)
                ranks_up += 1
            intruders = self.not_monophyletic.get(expected, []) if skipped else []
            self._writer.writerow(
                [
                    name,
                    expected,
                    taxon,
                    ranks_up,
                    new_clade or "",
                    ";".join(skipped),
                    len(intruders),
                    ";".join(intruders[:MAX_INTRUDERS]),
                ]
            )
            self.added += 1
            if skipped:
                entry = self.fallbacks.setdefault(expected, [0, taxon, ranks_up])
                entry[0] += 1
        self._file.flush()

    def summary(self) -> str:
        """Describe how many species were placed above their expected taxon, and which taxa caused it."""
        n_fallback = sum(count for count, _, _ in self.fallbacks.values())
        lines = []
        if not self.added:
            lines.append("No species were added.")
        elif not n_fallback:
            lines.append(f"Added {self.added} species, all to the lowest taxon with species in the backbone.")
        else:
            lines.append(
                f"Added {self.added} species; {n_fallback} ({n_fallback / self.added:.0%}) were placed in a "
                "higher taxon because the lowest taxon with species in the backbone is not monophyletic there."
            )
            ranked = sorted(self.fallbacks.items(), key=lambda item: (-item[1][0], item[0]))
            lines.append("Most affected:")
            for expected, (count, taxon, ranks_up) in ranked[:MAX_SUMMARY_TAXA]:
                intruders = self.not_monophyletic.get(expected, [])
                ranks = f"{ranks_up} rank{'s' if ranks_up != 1 else ''} up"
                broken_by = f"; broken up by {fmt_species_list(intruders)}" if intruders else ""
                lines.append(f"  {expected}: {count} species placed in {taxon} ({ranks}{broken_by})")
            if len(ranked) > MAX_SUMMARY_TAXA:
                lines.append(f"  and {len(ranked) - MAX_SUMMARY_TAXA} more taxa")
        lines.append(f"Placement of each added species: {self.path}")
        return "\n".join(lines)
