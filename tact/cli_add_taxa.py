#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Try to assign tips to a pre-existing tree based on a taxonomy
# Jonathan Chang, May 13, 2016

from __future__ import division
from __future__ import print_function

import csv
import logging
import operator
import random
import sys
from time import time

import click
import dendropy

from . import fastmrca
from .lib import crown_capture_probability
from .lib import get_new_times
from .tree_util import get_ages
from .tree_util import get_birth_death_rates
from .tree_util import get_min_age
from .tree_util import get_short_branches
from .tree_util import get_tip_labels
from .tree_util import graft_node
from .tree_util import is_binary
from .tree_util import is_fully_locked
from .tree_util import lock_clade
from .tree_util import update_tree_view
from .validation import validate_outgroups
from .validation import validate_taxonomy_tree
from .validation import BackboneCommand

logger = logging.getLogger(__name__)
# Speed up logging for PyPy
logging._srcfile = None
logging.logThreads = 0
logging.logProcesses = 0
logging.logMultiprocessing = 0

global invalid_map
invalid_map = {}

global mrca_rates
mrca_rates = {}


def search_ancestors_for_valid_backbone_node(taxonomy_node, backbone_tips, ccp):
    global invalid_map
    seen = []
    for anc in taxonomy_node.ancestor_iter():
        if anc.label in invalid_map:
            logger.debug(
                f"Cache HIT on invalid_map for {taxonomy_node.label} ({anc.label} => {invalid_map[anc.label].label})"
            )
            anc = invalid_map[anc.label]
        full_tax = get_tip_labels(anc)
        extant_tax = full_tax.intersection(backbone_tips)
        backbone_node = fastmrca.get(extant_tax)
        seen.append(anc.label)
        computed_ccp = crown_capture_probability(len(full_tax), len(extant_tax))
        if backbone_node is None:
            logger.info(f"    {taxonomy_node.label}: ancestor {anc.label} not monophyletic!")
        elif computed_ccp < ccp:
            logger.info(
                f"    {taxonomy_node.label}: using stem as ancestor {anc.label} ccp {computed_ccp:.2f} < {ccp}"
            )
            taxonomy_target = anc
            backbone_target = backbone_node.parent_node
            break
        else:
            taxonomy_target = anc
            backbone_target = backbone_node
            logger.info(f"    {taxonomy_node.label}: will instead assign these taxa to {taxonomy_target.label}")
            break
    else:
        anc_chain = " => ".join(seen)
        logger.error(f"Couldn't find valid taxonomy node in ancestor chain for {taxonomy_node.label} ({anc_chain})")
        sys.exit(1)
    seen.pop()  # ignore last node
    for x in seen:
        invalid_map[x] = taxonomy_target
    return (taxonomy_target, backbone_target)


def get_new_branching_times(backbone_node, taxonomy_node, told=None, tyoung=0, min_ccp=0.8, num_new_times=None):
    """
    Get `n_total` new branching times for a `node`.
    """
    global mrca_rates
    taxon = taxonomy_node.label
    birth, death, ccp, _ = mrca_rates[taxon]
    if ccp < min_ccp:
        if backbone_node.parent_node:
            new_told = backbone_node.parent_node.age
            if told is not None:
                logger.debug(f"    {taxon}: tmax {told:.2f} => {new_told:.2f} because ccp {ccp:.2f} < {min_ccp}")
            else:
                logger.debug(f"    {taxon}: tmax set to {new_told:.2f} because ccp {ccp:.2f} < {min_ccp}")
        else:
            # TODO: check for a root edge and graft a fake node above that
            new_told = backbone_node.age
            logger.debug(
                f"    {taxon}: tmax set to {new_told}; even though ccp {ccp:.2f} < {min_ccp}, clade is tree root"
            )
        told = new_told
    n_extant = len(backbone_node.leaf_nodes())
    n_total = len(taxonomy_node.leaf_nodes())
    if num_new_times is None:
        num_new_times = n_total - n_extant
    ages = get_ages(backbone_node)
    if len(backbone_node.leaf_nodes()) == 1 and told is None:
        # attach to stem in the case of a singleton
        told = backbone_node.parent_node.age
        logger.debug(f"    {taxon}: tmax set to {told:.2f} because taxon is singleton")
    if told is None:
        told = max(ages)
        logger.debug(f"    {taxon}: tmax set to {told:.2f} because of max age")
    logger.debug(
        f"    {taxon}: {num_new_times} new times: b={birth:.2f}, d={death:.2f}, tmax={told:.2f}, tmin={tyoung:.2f}"
    )
    times = get_new_times(ages, birth, death, num_new_times, told, tyoung)
    if len(times) > 5:
        logger.debug(f"    {taxon}: {times[0]:.2f}..{times[-1]:.2f}")
    else:
        logger.debug(f"    {taxon}: " + ", ".join(["{:.2f}" for x in times]).format(*times))
    return times


def fill_new_taxa(namespace, node, new_taxa, times, stem=False):
    for new_species, new_age in zip(new_taxa, times):
        new_node = dendropy.Node()
        new_node.annotations.add_new("creation_method", "fill_new_taxa")
        new_node.age = new_age
        new_leaf = new_node.new_child(taxon=namespace.require_taxon(new_species), edge_length=new_age)
        new_leaf.annotations.add_new("creation_method", "fill_new_taxa")
        new_leaf.age = 0
        node = graft_node(node, new_node, stem)

    count_short_branches = len(list(get_short_branches(node)))
    if count_short_branches:
        logger.info(f"{count_short_branches} short branches detected")

    return node


def create_clade(namespace, species, ages):
    tree = dendropy.Tree(taxon_namespace=namespace)
    species = list(species)
    ages.sort(reverse=True)
    # need to generate the "stem node"
    tree.seed_node.age = ages.pop(0)
    # clade of size 1?
    if not ages:
        node = tree.seed_node.new_child(edge_length=tree.seed_node.age, taxon=namespace.require_taxon(species[0]))
        node.age = 0.0
        for internal_node in tree.preorder_node_iter():
            internal_node.annotations.add_new("creation_method", "create_clade")
        return tree
    node = tree.seed_node.new_child()
    node.age = ages.pop(0)
    for age in ages:
        valid_nodes = [x for x in tree.nodes() if len(x.child_nodes()) < 2 and age < x.age and x != tree.seed_node]
        assert len(valid_nodes) > 0
        node = random.sample(valid_nodes, 1).pop()
        child = node.new_child()
        child.age = age
    n_species = len(species)
    random.shuffle(species)
    for node in tree.preorder_node_iter(filter_fn=lambda x: x.age > 0 and x != tree.seed_node):
        while len(node.child_nodes()) < 2 and len(species) > 0:
            new_species = species.pop()
            new_leaf = node.new_child(taxon=namespace.require_taxon(new_species))
            new_leaf.age = 0.0
    assert n_species == len(tree.leaf_nodes())
    assert len(tree.seed_node.child_nodes()) == 1
    for internal_node in tree.preorder_node_iter():
        internal_node.annotations.add_new("creation_method", "create_clade")
    assert is_binary(tree.seed_node.child_nodes()[0])
    tree.set_edge_lengths_from_node_ages(error_on_negative_edge_lengths=True)
    # Lock the child of the seed node so that things can still attach to the stem of this new clade
    lock_clade(tree.seed_node.child_nodes()[0])
    if list(get_short_branches(tree.seed_node)):
        logger.info("{} short branches detected".format(len(list(get_short_branches(tree.seed_node)))))
    return tree


def fmt_species_list(spp):
    spp = list(spp)
    if len(spp) > 2:
        return f"{spp[0]}, {spp[1]} and {len(spp) - 2} others"
    return " and ".join(spp)


def process_node(backbone_tree, backbone_bitmask, taxon_node, min_ccp, default_birth, default_death, yule=False):
    # TODO: Fix all the returns and refactor this into something sane
    global mrca_rates
    taxon = taxon_node.label
    parent = taxon_node.parent_node.label
    try:
        birth, death, ccp, _ = mrca_rates[parent]
    except KeyError:
        birth = default_birth
        death = default_death
        parent = "ROOT"
    species = get_tip_labels(taxon_node)
    if not taxon:
        logger.debug(f"MRCA: skipping unlabeled rank with {len(species)} species")
        return
    all_bitmask = backbone_tree.taxon_namespace.taxa_bitmask(labels=species)
    extant_bitmask = all_bitmask & backbone_bitmask
    if extant_bitmask is None or extant_bitmask == 0:
        logger.debug(f"MRCA: {taxon} not present in backbone")
        mrca_rates[taxon] = (birth, death, 0.0, f"from {parent} (unsampled)")
        return
    mrca = backbone_tree.mrca(leafset_bitmask=extant_bitmask)
    if not species.issuperset(get_tip_labels(mrca)):
        logger.debug(
            f"MRCA: {taxon} not monophyletic in backbone (from {fmt_species_list(get_tip_labels(mrca) - species)})"
        )
        mrca_rates[taxon] = (birth, death, 0.0, f"from {parent} (not monophyletic)")
        return
    extant = len(mrca.leaf_nodes())
    total = len(taxon_node.leaf_nodes())
    if extant > total:
        logger.warning(f"MRCA: {taxon} has {extant} extant species but should have {total} total species")
        mrca_rates[taxon] = (birth, death, 0, f"from {parent} (extant exceeds total)")
        return
    ccp = crown_capture_probability(total, extant)
    if total == 1:
        logger.debug(f"MRCA: {taxon} is a singleton")
        mrca_rates[taxon] = (birth, death, ccp, f"from {parent} (singleton)")
        return
    if total == 2:
        logger.debug(f"MRCA: {taxon} is a cherry")
        mrca_rates[taxon] = (birth, death, ccp, f"from {parent} (cherry)")
        return
    if ccp < min_ccp:
        logger.debug(
            f"MRCA: {taxon} has crown capture probability {ccp:.2f} < {min_ccp:.2f} ({extant}/{total} species)"
        )
        mrca_rates[taxon] = (birth, death, ccp, f"from {parent} (crown capture probability)")
        return
    sf = extant / total
    birth, death = get_birth_death_rates(mrca, sf, yule)
    logger.debug(f"MRCA: {taxon} b={birth:.2f}, d={death:.2f}, sf={sf:.2f} ({extant}/{total}), ccp={ccp:.2f}")
    mrca_rates[taxon] = (birth, death, ccp, "computed")


def run_precalcs(taxonomy_tree, backbone_tree, min_ccp=0.8, yule=False):
    global mrca_rates
    tree_tips = get_tip_labels(backbone_tree)
    backbone_bitmask = fastmrca.bitmask(tree_tips)
    all_possible_tips = get_tip_labels(taxonomy_tree)
    nnodes = len(taxonomy_tree.internal_nodes(exclude_seed_node=True))

    start_time = time()

    # Compute the rate of the root taxonomic node to use as a default value...
    logger.debug("Computing root birth and death rates.")
    extant_bitmask = backbone_bitmask & backbone_tree.taxon_namespace.taxa_bitmask(labels=all_possible_tips)
    root_mrca = backbone_tree.mrca(leafset_bitmask=extant_bitmask)
    root_birth, root_death = get_birth_death_rates(
        root_mrca, len(root_mrca.leaf_nodes()) / len(all_possible_tips), yule
    )

    with click.progressbar(
        taxonomy_tree.preorder_internal_node_iter(exclude_seed_node=True),
        width=12,
        label="Rates",
        length=nnodes,
        show_pos=True,
        item_show_func=lambda x: x.label if x else None,
    ) as progress:
        for node in progress:
            # updates global mrca_rates as a side effect
            process_node(backbone_tree, backbone_bitmask, node, min_ccp, root_birth, root_death, yule)

    diff = time() - start_time
    if diff > 1:
        logger.debug(f"FastMRCA calculation time: {diff:.1f} seconds")
    return mrca_rates


@click.command(cls=BackboneCommand)
@click.version_option(package_name="tact")
@click.option(
    "--taxonomy", help="a taxonomy tree", type=click.File("r"), required=True, callback=validate_taxonomy_tree
)
@click.option(
    "--backbone", help="the backbone tree to attach the taxonomy tree to", type=click.File("r"), required=True
)
@click.option(
    "--outgroups",
    help="comma separated list of outgroup taxa to ignore",
    type=click.UNPROCESSED,
    callback=validate_outgroups,
)
@click.option("--output", required=True, help="output base name to write out")
@click.option(
    "--min-ccp",
    default=0.8,
    type=click.FloatRange(0, 1, clamp=True),
    help="minimum probability to use to say that we've sampled the crown of a clade",
)
@click.option("--yule", help="assume a Yule pure-birth model (force extinction to be 0)", default=False, is_flag=True)
@click.option(
    "--ultrametricity-precision",
    default=1e-6,
    type=click.FloatRange(0, 1, clamp=True),
    help="precision for ultrametricity checks; by default, checks roughly digits of similarity",
)
@click.option("-v", "--verbose", help="emit extra information (can be repeated)", count=True)
def main(taxonomy, backbone, outgroups, output, min_ccp, verbose, yule, ultrametricity_precision):
    """
    Add tips onto a BACKBONE phylogeny using a TAXONOMY phylogeny.
    """
    logger.addHandler(logging.FileHandler(output + ".log.txt"))
    if verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
        logger.addHandler(logging.StreamHandler())

    tree = backbone
    tn = tree.taxon_namespace
    tree_tips = get_tip_labels(tree)
    all_possible_tips = get_tip_labels(taxonomy)

    logger.info(f"Backbone needs to add {len(tree_tips.symmetric_difference(all_possible_tips))} tips")

    full_clades = set()

    fastmrca.initialize(tree)

    with open(output + ".rates.csv", "w", encoding="utf-8") as wfile:
        writer = csv.writer(wfile)
        writer.writerow(("taxon", "birth", "death", "ccp", "source"))
        for key, value in run_precalcs(taxonomy, tree, min_ccp, yule=yule).items():
            row = [key]
            row.extend(value)
            writer.writerow(row)

    initial_length = len(tree_tips)

    bar = click.progressbar(
        label="TACT",
        length=len(all_possible_tips) - initial_length,
        show_pos=True,
        width=12,
        item_show_func=lambda x: x,
    )

    def bar_update():
        diff = len(tree_tips) - initial_length - bar.pos
        bar.update(diff, taxon if taxon else "")

    for taxon_node in taxonomy.postorder_internal_node_iter(exclude_seed_node=True):
        taxon = taxon_node.label
        if not taxon:
            continue
        species = get_tip_labels(taxon_node)
        extant_species = tree_tips.intersection(species)
        logger.info(f"**  {taxon} ({len(extant_species)}/{len(species)})  **")
        ccp = mrca_rates[taxon][2]

        clades_to_generate = full_clades.intersection(
            [x.label for x in taxon_node.postorder_internal_node_iter(exclude_seed_node=True)]
        )

        if not extant_species:
            # No species sampled, so create a clade from whole cloth
            logger.debug(f"    {taxon}: no species sampled, will create later")
            full_clades.add(taxon)
            bar_update()
            continue

        # Check for monophyly for this node
        node = fastmrca.get(extant_species)
        if not node:
            logger.info(f"    {taxon}: is not monophyletic")
            continue

        if extant_species == species:
            # Everything sampled and monophyletic, so skip this
            lock_clade(node)
            logger.debug(f"    {taxon}: all species accounted for")
            continue

        if tree_tips.issuperset(species):
            # XXX: Does this check ever get triggered?
            lock_clade(node)
            logger.info(f"    {taxon}: all species already present in tree")
            continue

        clade_ranks = [(clade, taxonomy.find_node_with_label(clade).level()) for clade in clades_to_generate]

        # Now add clades of unsampled species. Go from the lowest rank to
        # the highest (deepest level to lowest level). Shuffling before
        # sorting will randomize the order since Python uses stable sorting
        random.shuffle(clade_ranks)
        for clade, _ in sorted(clade_ranks, key=operator.itemgetter(1), reverse=True):
            full_node = taxonomy.find_node_with_label(clade)
            full_node_species = get_tip_labels(full_node)
            if tree_tips.issuperset(full_node_species):
                logger.info(f"    {taxon}: skipping clade {clade} as all species already present in tree")
                full_clades.remove(clade)
                continue
            logger.info(f"    {taxon}: adding clade {clade} (n={len(full_node.leaf_nodes())})")
            # Generate all times needed to attach to the main clade
            times = get_new_branching_times(
                node, taxon_node, tyoung=0, min_ccp=min_ccp, num_new_times=len(full_node_species)
            )

            if is_fully_locked(node):
                logger.info(f"    {taxon}: is fully locked, so attaching to stem")
                # Must attach to stem for this clade, so generate a time on the stem lineage
                times2 = get_new_branching_times(
                    node,
                    taxon_node,
                    min_ccp=min_ccp,
                    told=node.parent_node.age,
                    tyoung=node.age,
                    num_new_times=1,
                )
                # Drop the oldest time and add on our new time on the stem lineage
                times.sort()
                times.pop()
                times.append(times2.pop())
            else:
                # Even if the main clade isn't fully locked, it might have a constraint on a valid attachment point
                min_age = get_min_age(node)
                if min_age > 0 and max(times) < min_age:
                    logger.info(
                        f"    {taxon}: has a minimum age constraint {min_age:.2f}, "
                        + f"but oldest generated time was {max(times):.2f}"
                    )
                    times2 = get_new_branching_times(
                        node, taxon_node, tyoung=min_age, min_ccp=min_ccp, num_new_times=1
                    )
                    # Drop the oldest time and add on our new time on the stem lineage
                    times.sort()
                    times.pop()
                    times.append(times2.pop())

            # Generate a new tree
            new_tree = create_clade(tn, full_node_species, times)
            # Update our current MRCA node (because we might have attached to stem)
            node = graft_node(node, new_tree.seed_node, is_fully_locked(node) or ccp < min_ccp)
            tree_tips = update_tree_view(tree)
            # Update our view of what's in the tree
            extant_species = tree_tips.intersection(species)
            # We've added this clade so pop it off our stack
            full_clades.remove(clade)
            if not is_binary(node):
                raise ValueError("Tree is not binary!")
            bar_update()

        # Check to see if we need to continue adding species
        if extant_species == species:
            # Lock clade since it is monophyletic and filled
            lock_clade(node)
            # Skip taxon spray check
            continue
        if len(extant_species) == len(species):
            raise ValueError("Enough species are present but mismatched?")

        # Taxon spray
        logger.info(f"    {taxon}: adding {len(species.difference(extant_species))} new species")
        node = fastmrca.get(extant_species)
        times = get_new_branching_times(node, taxon_node, tyoung=get_min_age(node), min_ccp=min_ccp)
        node = fill_new_taxa(tn, node, species.difference(tree_tips), times, ccp < min_ccp)
        # Update stuff
        tree_tips = update_tree_view(tree)
        # Since only monophyletic nodes get to here, lock this clade
        lock_clade(node)
        if not is_binary(node):
            # Shouldn't happen
            raise ValueError("Tree is not binary!")
        bar_update()

    assert is_binary(tree.seed_node)
    # Reset terminal because we aren't using the context manager
    bar.render_finish()
    tree.ladderize()
    tree.write(path=output + ".newick.tre", schema="newick", suppress_rooting=True)
    tree.write(path=output + ".nexus.tre", schema="nexus")
    print()


if __name__ == "__main__":
    main()
