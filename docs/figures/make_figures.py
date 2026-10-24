"""Generate the SVG figures in docs/img.

Run from the repository root, in an environment where TACT is installed:

    poetry run python docs/figures/make_figures.py

The output is deterministic, so rerunning this without changing it leaves docs/img unchanged.
"""

from __future__ import annotations

import copy
import csv
import itertools
from collections.abc import Callable
from pathlib import Path

import dendropy
from svgfig import Figure, i, num
from trees import Node, draw_branches, from_newick, ladderize, layout, mark_unsampled, mrca, path_to

from tact.lib import crown_capture_probability

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "docs" / "img"

LEGEND_LOCKED = ("triangle", "clade", "ineligible for new grafts")
LEGEND_UNLOCKED = ("line", "unlocked", "eligible graft edges")
LEGEND_NEW = ("triangle", "new-clade", "newly generated clade")


def min_age_constraint() -> Figure:
    """Troubleshooting: a new clade whose attachment time is younger than every unlocked edge."""
    fig = Figure(
        660,
        380,
        "Minimum age constraint on a new clade",
        "A family whose sampled genera, Lycalopex and Vulpes, are complete and locked. Only their stems can "
        "receive grafts, and no part of those stems is younger than A, the crown age of Lycalopex. The unsampled "
        "genus Canis was given an attachment time B that is younger than A, where only locked edges exist, so TACT "
        "redraws the attachment time to be older than A and grafts Canis onto an unlocked stem.",
    )

    def x(age: float) -> float:
        return 10 + (11 - age) * 50

    crown, vulpes, a, b, graft, canis = 10, 7, 5, 2, 8, 1.5
    lyc = Node(a, label="Lycalopex")
    vul = Node(vulpes, label="Vulpes")
    lyc.cls = vul.cls = "unlocked"
    lyc.y, vul.y = 150, 250
    tree = Node(crown, [lyc, vul])
    layout(tree, x)

    fig.rect("band", x(crown), 56, x(a) - x(crown), 250)
    fig.text((x(crown) + x(a)) / 2, 296, "valid ages for attachment", "small ok-text", "middle")
    fig.path("ok-line dashed", f"M{num(x(a))},56 V306")
    fig.path("bad-line dashed", f"M{num(x(b))},56 V306")
    curve = f"C{num(x(b) - 50)},14 {num(x(graft) + 50)},14 {num(x(graft) + 1)},76"
    fig.arrow("bad-line", f"M{num(x(b) - 4)},52 {curve}")
    redraw = "B is too young, so redraw the attachment time ≥ A"
    fig.text((x(graft) + x(b)) / 2, 22, redraw, "small bad-text", "middle")

    draw_branches(fig, tree, "edge", stem_x=x(11))
    fig.circle("dot", tree.x, tree.y, 3.5)
    fig.text(tree.x, lyc.y - 12, "crown", "small muted", "middle")
    for leaf, half in ((lyc, 24), (vul, 28)):
        fig.triangle("clade", leaf.x, leaf.y, x(0), half)
        fig.padlock(x(0) - 24, leaf.y - 2)

    fig.path("new-edge", f"M{num(x(graft))},{num(lyc.y)} V84 H{num(x(canis))}")
    fig.triangle("new-clade", x(canis), 84, x(0), 14)
    fig.circle("new-dot", x(graft), lyc.y, 4.5)

    fig.text(x(0) + 12, 88, i("Canis"), "new-text")
    fig.text(x(0) + 12, lyc.y + 4, i("Lycalopex"))
    fig.text(x(0) + 12, vul.y + 4, i("Vulpes"))

    fig.path("axis", f"M{num(x(11))},322 H{num(x(0))} M{num(x(0))},316 V328")
    fig.path("ok-line", f"M{num(x(a))},316 V328")
    fig.path("bad-line", f"M{num(x(b))},316 V328")
    fig.text(x(a), 344, "A", "bold ok-text", "middle")
    fig.text(x(b), 344, "B", "bold bad-text", "middle")
    fig.text(x(0), 344, "present", "small muted", "middle")
    fig.text(x(11), 344, "older", "small muted")
    fig.legend(x(11), 366, [LEGEND_LOCKED, LEGEND_UNLOCKED, LEGEND_NEW])
    return fig


def fully_locked_stem() -> Figure:
    """Troubleshooting: a new clade grafted onto the stem of a clade whose every edge is locked."""
    fig = Figure(
        660,
        340,
        "Grafting onto the stem of a fully locked clade",
        "A subfamily Caninae in which only the genus Canis is sampled. Canis is complete and locked, so every edge "
        "below the sampled crown of Caninae is locked. TACT instead draws an attachment time on the stem of Canis, "
        "between its crown age and the age of its parent node, and grafts the unsampled genus Vulpes there. The "
        "graft point becomes the new crown of Caninae.",
    )

    def x(age: float) -> float:
        return 10 + (11 - age) * 50

    parent, crown, graft, vulpes = 10, 5, 7.5, 2.5
    canis = Node(crown, label="Canis")
    outgroup = Node(0, label="outgroup")
    canis.cls, outgroup.cls = "unlocked", "edge"
    canis.y, outgroup.y = 150, 260
    tree = Node(parent, [canis, outgroup])
    layout(tree, x)

    fig.rect("band", x(parent), 40, x(crown) - x(parent), 180)
    fig.text((x(parent) + x(crown)) / 2, 210, "valid ages on the stem", "small ok-text", "middle")

    draw_branches(fig, tree, "edge", stem_x=x(11))
    fig.circle("dot", tree.x, tree.y, 3.5)
    fig.text(tree.x + 10, tree.y - 9, "parent", "small muted")
    fig.text(x(0) + 12, outgroup.y + 4, "outgroup", "muted")

    fig.triangle("clade", canis.x, canis.y, x(0), 24)
    fig.padlock(x(0) - 24, canis.y - 2)
    fig.circle("dot", canis.x, canis.y, 3.5)
    fig.text((x(crown) + x(0)) / 2, 200, f"every edge in {i('Canis')} is ineligible", "small bad-text", "middle")

    fig.path("new-edge", f"M{num(x(graft))},{num(canis.y)} V84 H{num(x(vulpes))}")
    fig.triangle("new-clade", x(vulpes), 84, x(0), 14)
    fig.circle("new-dot", x(graft), canis.y, 4.5)
    fig.text(x(graft) + 9, 132, "new crown of Caninae", "small new-text")

    fig.text(x(0) + 12, 88, i("Vulpes"), "new-text")
    fig.text(x(0) + 12, canis.y + 4, i("Canis"))

    fig.path("axis", f"M{num(x(11))},280 H{num(x(0))} M{num(x(0))},275 V285")
    fig.path("ok-line", f"M{num(x(parent))},275 V285 M{num(x(crown))},275 V285")
    fig.text(x(parent), 300, "parent age", "small ok-text", "middle")
    fig.text(x(crown), 300, "crown age", "small ok-text", "middle")
    fig.text(x(0), 300, "present", "small muted", "middle")
    fig.legend(x(11), 320, [LEGEND_LOCKED, LEGEND_UNLOCKED, LEGEND_NEW])
    return fig


def intruder_taxon() -> Figure:
    """Troubleshooting: Cynoglossus inside the Soleidae clade of the tutorial backbone."""
    tree = dendropy.Tree.get(path=ROOT / "examples" / "Carangaria.tre", schema="newick", preserve_underscores=True)
    tree.is_rooted = True
    tree.calc_node_ages(is_force_max_age=True)
    with open(ROOT / "examples" / "Carangaria.csv", newline="") as fh:
        ranks = {row[-1].replace(" ", "_"): row for row in csv.reader(fh)}
    tips = [leaf.taxon.label for leaf in tree.leaf_node_iter()]

    def clade(column: int, name: str) -> dendropy.Node:
        return tree.mrca(taxon_labels=[tip for tip in tips if ranks[tip][column] == name])

    def has(node: dendropy.Node, column: int, name: str) -> bool:
        return any(ranks[leaf.taxon.label][column] == name for leaf in node.leaf_iter())

    pleuronectoidei = clade(2, "Pleuronectoidei")
    soleidae = clade(3, "Soleidae")
    parent = soleidae.parent_node
    (samaridae,) = (child for child in parent.child_node_iter() if child is not soleidae)
    (nested,) = (child for child in soleidae.child_node_iter() if has(child, 3, "Cynoglossidae"))
    (soles,) = (child for child in soleidae.child_node_iter() if child is not nested)
    (cynoglossus,) = (child for child in nested.child_node_iter() if has(child, 3, "Cynoglossidae"))
    (aseraggodes,) = (child for child in nested.child_node_iter() if child is not cynoglossus)
    first_step = next(child for child in pleuronectoidei.child_node_iter() if parent in child.preorder_iter())
    assert not has(soles, 3, "Cynoglossidae") and has(samaridae, 3, "Samaridae")
    assert {ranks[leaf.taxon.label][4] for leaf in aseraggodes.leaf_iter()} == {"Aseraggodes"}

    fig = Figure(
        800,
        400,
        "An intruding taxon breaks up a family",
        "Part of the Carangaria backbone from the tutorial. The smallest clade containing every sampled species of "
        "Soleidae also contains Cynoglossus, from the family Cynoglossidae, so Soleidae is not monophyletic. TACT "
        "therefore adds the unsampled Soleidae genus Brachirus anywhere in the enclosing suborder Pleuronectoidei, "
        "for example next to Samaridae, rather than within Soleidae.",
        "Node ages come from examples/Carangaria.tre; the nodes between the Pleuronectoidei crown and the "
        "Samaridae + Soleidae clade are collapsed, and the position of Brachirus is one example.",
    )

    def x(age: float) -> float:
        return 60 + (75 - age) * 500 / 75

    def leaf(node: dendropy.Node, y: float) -> Node:
        tip = Node(node.age)
        tip.y = y
        return tip

    sole_a, cyn, ase = leaf(soles, 60), leaf(cynoglossus, 120), leaf(aseraggodes, 170)
    sam, rest = leaf(samaridae, 225), leaf(first_step, 325)
    s = Node(soleidae.age, [sole_a, Node(nested.age, [cyn, ase])])
    top = Node(pleuronectoidei.age, [Node(parent.age, [s, sam]), rest])
    layout(top, x)
    graft, brachirus = 40, 8  # an example placement

    fig.rect("sole-box", s.x - 8, 34, 560, 152, 6)
    fig.text(s.x + 2, 28, "smallest clade containing every sampled species of Soleidae", "small sole-text")
    draw_branches(fig, top, "unlocked", stem_x=x(75) - 10)
    for node, cls, half in ((sole_a, "sole", 16), (cyn, "intruder", 14), (ase, "sole", 7)):
        fig.triangle(cls, node.x, node.y, x(0), half)
    fig.triangle("clade", sam.x, sam.y, x(0), 10)
    fig.triangle("clade", rest.x, rest.y, x(0), 25)

    fig.path("new-edge", f"M{num(x(graft))},{num(sam.y)} V268 H{num(x(brachirus))}")
    fig.triangle("new-clade", x(brachirus), 268, x(0), 9)
    fig.circle("new-dot", x(graft), sam.y, 4.5)

    lx = x(0) + 12
    fig.text(lx, sole_a.y - 2, f"{i('Solea')}, {i('Pegusa')} and others", "sole-text")
    fig.text(lx, sole_a.y + 14, "Soleidae", "small sole-text")
    fig.text(lx, cyn.y - 2, i("Cynoglossus"), "bad-text")
    fig.text(lx, cyn.y + 14, "Cynoglossidae", "small bad-text")
    fig.text(lx, ase.y + 4, i("Aseraggodes"), "sole-text")
    fig.text(lx, sam.y + 4, "Samaridae", "muted")
    fig.text(lx, 272, f'{i("Brachirus")} <tspan class="small">(new, Soleidae)</tspan>', "new-text")
    fig.text(lx, rest.y + 4, "other Pleuronectoidei", "muted")

    mid = (sole_a.y - 16 + rest.y + 25) / 2
    fig.path("bracket", f"M760,{num(sole_a.y - 16)} h6 V{num(rest.y + 25)} h-6 M766,{num(mid)} h6")
    fig.add(
        f'<text x="780" y="{num(mid)}" class="small muted" transform="rotate(90 780 {num(mid)})" '
        'text-anchor="middle">Pleuronectoidei</text>'
    )
    fig.legend(50, 382, [("line", "unlocked", f"where TACT can place {i('Brachirus')}"), LEGEND_NEW,
                         ("triangle", "intruder", "intruding taxon")])  # fmt: skip
    return fig


# Newick trees from fig1.R, which drew Figure 1 of Chang et al. (2020)
COMPLETE_TREE = (
    "((((((c10:0.06293181943,c9:0.06293181943):1.087444137,(c8:0.2746875215,c7:0.2746875215):0.875688435)"
    ":0.1630523222,c6:1.313428279):1.226515422,((c5:0.4823147685,c4:0.4823147685):0.6806888812,c3:1.16300365)"
    ":1.376940051):1.108085182,(c2:0.9160338094,c1:0.9160338094):2.731995073):1.342196861,((((((b6:0.3074985473,"
    "b5:0.3074985473):0.1467284868,b4:0.454227034):1.089381333,b3:1.543608367):0.4342465053,b2:1.977854872)"
    ":1.528606819,b1:3.506461691):0.8984834813,(((a4:2.217566675,a3:2.217566675):1.018220869,a2:3.235787545)"
    ":1.007163139,a1:4.242950684):0.1619944884):0.5852805708);"
)
TACT_TREE = (
    "((((c10:1.150375956,(((c8:0.09628559172,c7:0.09628559252):0.5074804988,c6:0.6037660916):0.2445160628,"
    "c9:0.8482821544):0.3020938034):1.389567744,(c5:1.346628605,c4:1.346628605):1.193315096):1.108085182,"
    "((c1:0.4772114505,c2:0.477211452):0.3805352588,c3:0.8577467108):2.790282173):1.342196861,((((b5:0.4542270341,"
    "b4:0.454227034):1.089381333,b3:1.543608367):1.962853324,((b1:0.9099970704,b2:0.909997071):0.7500736044,"
    "b6:1.660070675):1.846391017):0.8984834813,(((a4:2.217566675,a3:2.217566675):1.018220869,a2:3.235787545)"
    ":1.007163139,a1:4.242950684):0.1619944884):0.5852805708);"
)
UNSAMPLED = {"b2", "b6", "c2", "c3", "c4", "c7", "c6", "c9"}


def incomplete_sampling() -> Figure:
    """Background: the complete tree, sampling fraction, taxonomic and TACT views of one tree."""
    panel_w, x0, tree_w, y0, dy = 200, 24, 128, 58, 10.5
    complete = from_newick(COMPLETE_TREE)
    ladderize(complete)
    tmax = complete.age
    height = y0 + len(complete.leaves()) * dy + 40
    fig = Figure(
        4 * panel_w + 20,
        height,
        "Ways to handle an incompletely sampled phylogeny",
        "Four versions of the same tree of 20 species in three clades of 10, 6 and 4 species. (a) The complete "
        "tree. (b) The sampled tree, with the unsampled lineages dashed and the fraction of each clade sampled. "
        "(c) The taxonomic method, with each clade collapsed to a single tip with its species richness. (d) One "
        "tree from TACT, with the unsampled species placed at random within their clades.",
        "Trees from fig1.R, which drew Figure 1 of Chang et al. (2020).",
    )

    def panel(index: int, letter: str, heading: str) -> tuple[float, float]:
        left = x0 + index * panel_w
        fig.text(left, 30, f'<tspan class="bold">{letter}</tspan>  {heading}')
        return left, left + 10 + tree_w

    def scale(left: float) -> Callable[[float], float]:
        return lambda age: left + 10 + (tmax - age) * tree_w / tmax

    def clade_tips(tree: Node, prefix: str) -> list[Node]:
        return [leaf for leaf in tree.leaves() if leaf.label.startswith(prefix)]

    # a) complete tree
    left, _ = panel(0, "a", "complete tree")
    layout(complete, scale(left), y0, dy)
    draw_branches(fig, complete, stem_x=left)

    # b) sampling fraction: only the sampled lineages, with the fraction sampled in each clade
    left, right = panel(1, "b", "sampling fraction")
    sampled = copy.deepcopy(complete)
    for leaf in sampled.leaves():
        leaf.sampled = leaf.label not in UNSAMPLED
    mark_unsampled(sampled, "unsampled")
    layout(sampled, scale(left), y0, dy)
    draw_branches(fig, sampled, stem_x=left)
    for prefix in "cba":
        tips = clade_tips(sampled, prefix)
        lo, hi = tips[0].y, tips[-1].y
        fraction = sum(tip.sampled for tip in tips) / len(tips)
        fig.path("bracket", f"M{num(right + 8)},{num(lo - 3)} V{num(hi + 3)}")
        fig.text(right + 14, (lo + hi) / 2 + 4, f'<tspan class="italic">f</tspan> = {fraction:.2g}', "small fg")

    # c) taxonomic: each clade collapsed to a tip with its number of species
    left, right = panel(2, "c", "taxonomic")
    leaves = {}
    for prefix in "cba":
        tips = clade_tips(complete, prefix)
        node = Node(mrca(complete, {tip.label for tip in tips}).age, label=prefix)
        node.y = (tips[0].y + tips[-1].y) / 2
        leaves[prefix] = (node, (tips[-1].y - tips[0].y) / 2 + 3, len(tips))
    ba = mrca(complete, {tip.label for tip in clade_tips(complete, "a") + clade_tips(complete, "b")})
    skeleton = Node(tmax, [leaves["c"][0], Node(ba.age, [leaves["b"][0], leaves["a"][0]])])
    layout(skeleton, scale(left))
    draw_branches(fig, skeleton, stem_x=left)
    for node, half, total in leaves.values():
        fig.triangle("clade", node.x, node.y, right, half)
        fig.text(right + 8, node.y + 4, f'<tspan class="italic">N</tspan> = {total}', "small fg")

    # d) TACT: the unsampled species placed at random within their clades
    left, _ = panel(3, "d", "TACT")
    tact = from_newick(TACT_TREE)
    ladderize(tact)
    for leaf in tact.leaves():
        leaf.sampled = leaf.label not in UNSAMPLED
    mark_unsampled(tact, "placed")
    layout(tact, scale(left), y0, dy)
    draw_branches(fig, tact, stem_x=left)

    fig.legend(
        x0,
        height - 16,
        [
            ("line", "branch", "sampled lineage"),
            ("line", "unsampled", "unsampled lineage"),
            ("line", "placed", "species placed by TACT"),
        ],
    )
    return fig


# A genus of 10 species whose crown splits it into clades of 6 and 4
GENUS_TREE = (
    "((t1:4.500,(t2:3.563,(t3:0.868,t4:0.868):2.695):0.937):3.000,"
    "((t5:3.956,(t6:0.985,t7:0.985):2.971):1.544,(t8:2.622,(t9:2.194,t10:2.194):0.428):2.878):2.000);"
)


def crown_capture() -> Figure:
    """Troubleshooting: a sample of 3 species that does, and does not, include the crown of a genus of 10."""
    panel_w, y0, dy, stem_age = 400, 70, 22, 10.0
    genus = from_newick(GENUS_TREE)
    ladderize(genus, larger_first=False)
    n, k = len(genus.leaves()), 3
    captured = crown_capture_probability(n, k)
    height = y0 + (n - 1) * dy + 80
    fig = Figure(
        2 * panel_w,
        height,
        "Crown capture probability",
        f"The same genus of {n} species, with {k} species sampled in two different ways. (a) The sampled species "
        "include at least one on each side of the true crown, so the crown of the sampled tree is the true crown. "
        "(b) All the sampled species are on one side of the true crown, so the sampled crown is younger than the "
        "true crown, and the unsampled species on the other side join the tree on the stem of the sampled crown.",
        f"Probabilities from tact.lib.crown_capture_probability({n}, {k}).",
    )

    def panel(index: int, picks: set[str], letter: str, heading: str, probability: float) -> None:
        left = index * panel_w
        tree = copy.deepcopy(genus)
        for leaf in tree.leaves():
            leaf.sampled = leaf.label in picks
        mark_unsampled(tree, "unsampled")
        x0 = left + 30
        layout(tree, lambda age: x0 + (stem_age - age) * 30, y0, dy)
        fig.text(left + 20, 30, f'<tspan class="bold">{letter}</tspan>  {heading}')
        draw_branches(fig, tree, stem_x=x0)
        sampled_crown = mrca(tree, picks)
        if sampled_crown is not tree:
            # the stem of the sampled crown runs back to where the genus joins the rest of the tree
            route = path_to(tree, sampled_crown)
            d = f"M{num(x0)},{num(tree.y)} H{num(tree.x)}"
            d += "".join(f" M{num(a.x)},{num(a.y)} V{num(b.y)} H{num(b.x)}" for a, b in itertools.pairwise(route))
            fig.path("unlocked", d)
            fig.text(tree.x - 10, tree.y - 8, "stem of the sampled crown", "small ok-text", "end")
        for leaf in tree.leaves():
            if leaf.sampled:
                fig.circle("tip", leaf.x, leaf.y, 3.5)
        if sampled_crown is tree:
            fig.circle("crown", tree.x, tree.y, 5)
            fig.text(
                tree.x - 10, tree.y + 20, f'sampled crown<tspan x="{num(tree.x - 10)}" dy="15">= true crown</tspan>',
                "small fg", "end",
            )  # fmt: skip
        else:
            fig.circle("true-crown", tree.x, tree.y, 5)
            fig.text(tree.x - 10, tree.y + 20, "true crown", "small muted", "end")
            fig.circle("crown", sampled_crown.x, sampled_crown.y, 5)
            fig.text(sampled_crown.x - 10, sampled_crown.y + 20, "sampled crown", "small fg", "end")
        fig.text(left + 20, y0 + (n - 1) * dy + 38, f"probability {probability:.2f} with {n} species, {k} sampled",
                 "small muted")  # fmt: skip

    panel(0, {"t2", "t6", "t9"}, "a", "the sample spans the true crown", captured)
    panel(1, {"t6", "t8", "t10"}, "b", "the sample misses the true crown", 1 - captured)
    fig.legend(
        20,
        height - 16,
        [
            ("tip", "branch", "sampled species"),
            ("line", "unsampled", "unsampled species"),
            ("line", "unlocked", "stem of the sampled crown"),
        ],
    )
    return fig


FIGURES = {
    "min-age-constraint": min_age_constraint,
    "fully-locked-stem": fully_locked_stem,
    "intruder-taxon": intruder_taxon,
    "crown-capture": crown_capture,
    "incomplete-sampling": incomplete_sampling,
}


def main() -> None:
    """Write every figure to docs/img."""
    for name, make in FIGURES.items():
        path = OUT / f"{name}.svg"
        make().write(path)
        print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
