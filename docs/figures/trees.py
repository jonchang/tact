"""Rectangular tree layout for the documentation figures."""

from __future__ import annotations

import math
from collections.abc import Callable, Iterator

import dendropy
from svgfig import Figure, num


class Node:
    """A node in a tree to draw, with its age and the class of the branch leading to it."""

    def __init__(self, age: float = 0.0, children: list[Node] | None = None, label: str = "") -> None:
        """Make a node; a node without children is a tip."""
        self.age = age
        self.children = children or []
        self.label = label
        self.cls: str | None = None
        self.sampled = True
        self.x = 0.0
        self.y = math.nan  # set by layout

    def leaves(self) -> list[Node]:
        """The tips below this node, in drawing order."""
        if not self.children:
            return [self]
        return [leaf for child in self.children for leaf in child.leaves()]

    def nodes(self) -> Iterator[Node]:
        """This node and every node below it, in preorder."""
        yield self
        for child in self.children:
            yield from child.nodes()

    def find(self, label: str) -> Node:
        """The node with this label."""
        return next(node for node in self.nodes() if node.label == label)


def from_newick(newick: str) -> Node:
    """Read a tree with branch lengths, setting each node's age from the tips."""
    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    tree.calc_node_ages(is_force_max_age=True)

    def convert(node: dendropy.Node) -> Node:
        label = node.taxon.label if node.taxon else (node.label or "")
        return Node(node.age, [convert(child) for child in node.child_node_iter()], label)

    return convert(tree.seed_node)


def ladderize(node: Node, larger_first: bool = True) -> None:
    """Order children by size, larger clades first (drawn above) unless `larger_first`; ties keep their order."""
    for child in node.children:
        ladderize(child, larger_first)
    node.children.sort(key=lambda child: len(child.leaves()) * (-1 if larger_first else 1))


def mrca(root: Node, labels: set[str]) -> Node:
    """The most recent common ancestor of the tips with these labels."""
    node = root
    while True:
        for child in node.children:
            if labels <= {leaf.label for leaf in child.leaves()}:
                node = child
                break
        else:
            return node


def path_to(root: Node, target: Node) -> list[Node]:
    """The nodes from `root` down to `target`, inclusive."""
    if root is target:
        return [root]
    for child in root.children:
        path = path_to(child, target)
        if path:
            return [root, *path]
    return []


def mark_unsampled(node: Node, cls: str) -> bool:
    """Give the class `cls` to branches with no sampled tip below them; returns whether any tip below is sampled."""
    if node.children:
        sampled = False
        for child in node.children:
            sampled |= mark_unsampled(child, cls)
    else:
        sampled = node.sampled
    if not sampled:
        node.cls = cls
    return sampled


def layout(root: Node, x_of_age: Callable[[float], float], y0: float = 0.0, dy: float = 10.0) -> None:
    """Place nodes: x from age, tips dy apart from y0 unless already placed, and parents midway between children."""
    y = y0
    for leaf in root.leaves():
        if math.isnan(leaf.y):
            leaf.y = y
        y = leaf.y + dy

    def place(node: Node) -> None:
        node.x = x_of_age(node.age)
        if node.children:
            for child in node.children:
                place(child)
            node.y = (node.children[0].y + node.children[-1].y) / 2

    place(root)


def draw_branches(fig: Figure, root: Node, default: str = "branch", stem_x: float | None = None) -> None:
    """Draw every branch as a corner from its parent, one path per class."""
    by_cls: dict[str, list[str]] = {}
    if stem_x is not None:
        by_cls.setdefault(root.cls or default, []).append(f"M{num(stem_x)},{num(root.y)} H{num(root.x)}")
    for node in root.nodes():
        for child in node.children:
            d = f"M{num(node.x)},{num(node.y)} V{num(child.y)} H{num(child.x)}"
            by_cls.setdefault(child.cls or default, []).append(d)
    for cls, parts in by_cls.items():
        fig.path(cls, " ".join(parts))
