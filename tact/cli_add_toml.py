"""Command line interface to assign tips to a pre-existing tree based on a TOML configuration file."""

from __future__ import annotations

import copy
import logging
import operator
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from functools import reduce
from typing import IO, Any

import click
import dendropy
import toml

from .console import progressbar, setup_logging
from .lib import get_new_times
from .tree_util import (
    get_ages,
    get_birth_death_rates,
    get_min_age,
    get_tip_labels,
    graft_node,
    is_binary,
    lock_clade,
    node_age,
    unlock_clade,
    update_tree_view,
)
from .validation import BackboneCommand

logger = logging.getLogger(__name__)


@dataclass
class TactConstraint:
    """Class for keeping track of a constraint in TACT (positive or negative)."""

    mrca: list[str] = field(default_factory=list)
    stem: bool = False

    def __post_init__(self) -> None:
        self.mrca = [x.replace("_", " ") for x in self.mrca]


@dataclass(repr=False)
class TactItem:
    """Class for keeping track of things to TACT."""

    name: str
    missing: int
    include: list[TactConstraint] = field(default_factory=list)
    exclude: list[TactConstraint] = field(default_factory=list)
    preserve_generic_monophyly: bool = False

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> TactItem:
        """Build a TactItem from a `[[tact]]` table in the TOML configuration."""
        config = dict(config)
        config["include"] = [TactConstraint(**x) for x in config.get("include") or []]
        config["exclude"] = [TactConstraint(**x) for x in config.get("exclude") or []]
        return cls(**config)

    def __repr__(self) -> str:
        return (
            f"TactItem('{self.name}', missing={self.missing}, "
            + f"include={len(self.include)}, exclude={len(self.exclude)})"
        )

    def __post_init__(self) -> None:
        if not self.include:
            logger.error(f"Need at least one include for {self.name}!")
            sys.exit(1)

        # Check restriction where singletons with include must have stem = True
        for constraint in self.include:
            if len(constraint.mrca) == 1 and constraint.stem is False:
                logger.error(f"Include specifications for singleton invalid without `stem = true`:\n{constraint}")
                sys.exit(1)


def ensure_mrca(tree: dendropy.Tree, tips: list[str], node: dendropy.Node | None = None) -> dendropy.Node:
    """Perform initial checks to ensure we can do MRCA calculations."""
    try:
        node = node if node else tree.seed_node
        return tree.mrca(taxon_labels=tips, start_node=node)
    except KeyError:
        diff = set(tips) - set(get_tip_labels(tree))
        if diff:
            logger.error("Could not find MRCA, these tips were not present in tree:")
            logger.error(list(diff))
        else:
            logger.error("Could not find MRCA for these tips:")
            logger.error(tips)
        logger.error("Ensure that:")
        logger.error("* your tips are properly spelled and present in the tree")
        if node != tree.seed_node:
            logger.error("* they are properly nested in `include` specifications")
        sys.exit(1)


def do_tact(tree: dendropy.Tree, item: TactItem) -> dendropy.Tree:
    # First, get the MRCA of _all_ `include` leafs. This is the basis of our rate computation,
    # and how we actually implement polyphyletic groups.
    included_tips = reduce(operator.iadd, [x.mrca for x in item.include], [])
    mrca_node = ensure_mrca(tree, included_tips)

    # Compute the rates on that (possibly expansive) MRCA node.
    extant_tips = len(mrca_node.leaf_nodes())

    should_include_root = extant_tips == 1 and len(item.include) == 1

    birth, death = get_birth_death_rates(
        mrca_node, extant_tips / (extant_tips + item.missing), include_root=should_include_root
    )
    logger.info(f"{item.name} => b={birth}, d={death}")

    # First, lock everything descending from the MRCA node, including its stem
    lock_clade(mrca_node, True)

    # Then, selectively unlock `include`s.
    for include in item.include:
        # If this was a single include with stem=True, then it gets unlocked here.
        inner_mrca_node = ensure_mrca(tree, include.mrca, mrca_node)
        unlock_clade(inner_mrca_node, include.stem)
        # Generic monophyly is ensured by pretending these are `exclude` clades,
        # except that we don't exclude them if the leaves of the current include
        # node belong wholly to that genus.
        if item.preserve_generic_monophyly:
            genera_map = defaultdict(set)
            for tip in get_tip_labels(inner_mrca_node):
                genus, _ = re.split("[_ ]+", tip, maxsplit=1)
                genera_map[genus].add(tip)

            if len(genera_map) > 1:
                for species in genera_map.values():
                    node = tree.mrca(taxon_labels=species, start_node=inner_mrca_node)
                    if node and species == get_tip_labels(node):
                        if len(species) == 1:
                            # Lock the stem of singletons
                            lock_clade(node, stem=True)
                        else:
                            lock_clade(node)

    # Finally, re-lock `exclude`s
    for exclude in item.exclude:
        lock_clade(ensure_mrca(tree, exclude.mrca, mrca_node), exclude.stem)

    # Compute new branching times
    ages = get_ages(mrca_node, include_root=should_include_root)
    times = get_new_times(ages, birth, death, item.missing, max(ages), get_min_age(mrca_node))

    # TODO: currently does not account for the possibility of a disjoint set of edges.
    # If this happens, reroll the time?
    for idx, time in enumerate(times):
        new_name = f"{item.name} tact {idx}"
        logger.info(f"=> Grafting {new_name} @ {time}")
        new_node = dendropy.Node()
        new_node.age = time
        new_leaf = new_node.new_child(taxon=dendropy.Taxon(label=new_name), edge_length=time)
        new_leaf.age = 0
        # Assume stem is fair game, since it would have been unlocked or kept locked earlier
        mrca_node = graft_node(mrca_node, new_node, True)

    tree.reconstruct_taxon_namespace()
    update_tree_view(tree)
    return tree


def do_replicate(backbone: dendropy.Tree, to_tact: list[TactItem], label: int) -> dendropy.Tree:
    """Perform a replicate of a TACT analysis."""
    logger.info(f"<<< Replicate {label} >>>")
    tree = copy.deepcopy(backbone)
    for item in to_tact:
        tree = do_tact(tree, item)
    if not is_binary(tree.seed_node):
        logger.error("Tree is not binary!")
        sys.exit(1)
    tree.ladderize()
    return tree


@click.command(cls=BackboneCommand)
@click.option("--config", help="configuration file", type=click.File("r"), required=True)
@click.option("--backbone", help="the backbone tree", type=click.File("r"), required=True)
@click.option("--output", required=True, help="output base name to write out")
@click.option("-v", "--verbose", help="emit extra information (can be repeated)", count=True)
@click.option(
    "--ultrametricity-precision",
    help="precision for ultrametricity checks; by default, checks roughly digits of similarity",
    default=1e-6,
    type=click.FloatRange(0, 1, clamp=True),
)
@click.option("--replicates", help="how many tacted trees to create", default=10, type=click.IntRange(1))
@click.option(
    "--cores",
    help="how many parallel cores to use",
    type=click.IntRange(1, os.cpu_count() or 1, clamp=True),
    default=os.cpu_count() or 1,
)
def main(
    config: IO[str],
    backbone: dendropy.Tree,
    output: str,
    verbose: int,
    ultrametricity_precision: float,
    replicates: int,
    cores: int,
) -> None:
    """Add tips onto a BACKBONE phylogeny using a CONFIG file."""
    setup_logging(output, verbose)

    logger.info("Reading configuration")
    settings = toml.load(config)

    to_tact = [TactItem.from_config(x) for x in settings["tact"]]

    # Ensure the proper ordering of TACT items based on divergence time of implied MRCA nodes
    to_tact.sort(
        key=lambda ii: node_age(ensure_mrca(backbone, reduce(operator.iadd, [x.mrca for x in ii.include], [])))
    )

    # Compute global birth/death rates. Not currently used (but could be?)
    backbone_tips = len(backbone.leaf_nodes())
    root_sf = backbone_tips / (backbone_tips + sum([x.missing for x in to_tact]))
    root_birth, root_death = get_birth_death_rates(backbone.seed_node, root_sf)
    logger.info(f"Root b={root_birth}, d={root_death}, sf={root_sf}")

    results = []
    with ProcessPoolExecutor(max_workers=cores) as executor:
        acc = []
        for idx in range(replicates):
            acc.append(executor.submit(do_replicate, backbone, to_tact, idx))
        with progressbar(as_completed(acc), label="TACT", length=replicates) as rf:
            for future in rf:
                results.append(future.result())

    ntip = len(results[0].leaf_nodes())
    for result_tree in results:
        new_ntip = len(result_tree.leaf_nodes())
        if new_ntip != ntip:
            logger.warning("TACTed trees have differing numbers of tips!")
            logger.warning(f"{ntip} != {new_ntip}")

    forest = dendropy.TreeList(results)
    forest.write(path=output + ".newick.tre", schema="newick", suppress_rooting=True)
    forest.write(path=output + ".nexus.tre", schema="nexus")


if __name__ == "__main__":
    main()
