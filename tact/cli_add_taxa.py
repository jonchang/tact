#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Try to assign tips to a pre-existing tree based on a taxonomy
# Jonathan Chang, May 13, 2016

# Pragmata
from __future__ import division
from __future__ import print_function

# Internal
from .lib import get_birth_death_rates, get_ages, is_binary, get_short_branches, get_tip_labels, crown_capture_probability, edge_iter, get_new_times
from . import fastmrca

# Python standard library
import logging
import operator
import random
import sys
from time import time
import collections
import csv

# Third party
import dendropy
import click

logger = logging.getLogger(__name__)
# Speed up logging for pypy
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
            logger.debug("Cache HIT on invalid_map for {} ({} => {})".format(taxonomy_node.label, anc.label, invalid_map[anc.label].label))
            anc = invalid_map[anc.label]
        full_tax = get_tip_labels(anc)
        extant_tax = full_tax.intersection(backbone_tips)
        backbone_node = fastmrca.get(extant_tax)
        seen.append(anc.label)
        if backbone_node is None:
            logger.info("    {}: ancestor {} not monophyletic!".format(taxonomy_node.label, anc.label))
        elif crown_capture_probability(len(full_tax), len(extant_tax)) < ccp:
            logger.info("    {}: ancestor {} fails crown threshold ({} < {}); using stem".format(taxonomy_node.label, anc.label, crown_capture_probability(len(full_tax), len(extant_tax)), ccp))
            taxonomy_target = anc
            backbone_target = backbone_node.parent_node
            break
        else:
            taxonomy_target = anc
            backbone_target = backbone_node
            logger.info("    {}: will instead assign these taxa to {}".format(taxonomy_node.label, taxonomy_target.label))
            break
    else:
        logger.error("Couldn't find valid taxonomy node in ancestor chain for {} ({})".format(taxonomy_node.label, " => ".join(seen)))
        sys.exit(1)
    seen.pop()  # ignore last node
    for x in seen:
        invalid_map[x] = taxonomy_target
    return (taxonomy_target, backbone_target)


def get_new_branching_times(backbone_node, taxonomy_node, backbone_tree, told=None, tyoung=0, min_ccp=0.8, num_new_times=None):
    """
    Get `n_total` new branching times for a `node`.
    """
    global mrca_rates
    taxon = taxonomy_node.label
    birth, death, ccp, source = mrca_rates[taxon]
    if ccp < min_ccp:
        if backbone_node.parent_node:
            new_told = backbone_node.parent_node.age
            if told is not None:
                logger.debug("    {}: tmax {} => {} because ccp {:.2f} < {}".format(taxon, told, new_told, ccp, min_ccp))
            else:
                logger.debug("    {}: tmax set to {} because ccp {:.2f} < {}".format(taxon, new_told, ccp, min_ccp))
        else:
            # TODO: check for a root edge and graft a fake node above that
            new_told = backbone_node.age
            logger.debug("    {}: tmax set to {} because even though ccp {:.2f} < {} clade is tree root".format(taxon, new_told, ccp, min_ccp))
        told = new_told
    n_extant = len(backbone_node.leaf_nodes())
    n_total = len(taxonomy_node.leaf_nodes())
    if num_new_times is None:
        num_new_times = n_total - n_extant
    ages = get_ages(backbone_node)
    if len(backbone_node.leaf_nodes()) == 1 and told is None:
        # attach to stem in the case of a singleton
        told = backbone_node.parent_node.age
        logger.debug("    {}: tmax set to {:.2f} because taxon is singleton".format(taxon, told))
    if told is None:
        told = max(ages)
        logger.debug("    {}: tmax set to {:.2f} because of max age".format(taxon, told))
    logger.debug("    {}: {} new times: b={:.2f}, d={:.2f}, tmax={:.2f}, tmin={:.2f}".format(taxon, num_new_times, birth, death, told, tyoung))
    times = get_new_times(ages, birth, death, num_new_times, told, tyoung)
    if len(times) > 5:
        logger.debug("    {}: {:.2f}..{:.2f}".format(taxon, times[0], times[-1]))
    else:
        logger.debug(("    {}: " + ", ".join(["{:.2f}" for x in times])).format(taxon, *times))
    return times


def fill_new_taxa(namespace, node, new_taxa, times, stem=False, excluded_nodes=None):
    for new_species, new_age in zip(new_taxa, times):
        new_node = dendropy.Node()
        new_node.annotations.add_new("creation_method", "fill_new_taxa")
        new_node.age = new_age
        new_leaf = new_node.new_child(taxon=namespace.require_taxon(new_species), edge_length=new_age)
        new_leaf.annotations.add_new("creation_method", "fill_new_taxa")
        new_leaf.age = 0
        node = graft_node(node, new_node, stem)

    if list(get_short_branches(node)):
        logger.warning("{} short branches detected".format(len(list(get_short_branches(node)))))

    return node


def graft_node(graft_recipient, graft, stem=False):
    """
    Grafts a node `graft` randomly in the subtree below node
    `graft_recipient`. The attribute `graft.age` must be set so
    we know where is the best place to graft the node. The node
    `graft` can optionally have child nodes, in this case the
    `edge.length` attribute should be set on all child nodes if
    the tree is to remain ultrametric.
    """

    # We graft things "below" a node by picking one of the children
    # of that node and forcing it to be sister to the grafted node
    # and adjusting the edge lengths accordingly. Therefore, the node
    # *above* which the graft lives (i.e., the one that will be the child
    # of the new graft) must fulfill the following requirements:
    #
    # 1. Must not be the crown node (cannot graft things above crown node)
    # 2. Must be younger than the graft node (no negative branches)
    # 3. Seed node must be older than graft node (no negative branches)
    # 4. Must not be locked (intruding on monophyly)
    def filter_fn(x):
        return x.head_node.age <= graft.age and x.head_node.parent_node.age >= graft.age and x.label != "locked"
    all_edges = list(edge_iter(graft_recipient))
    if stem:
        # also include the crown node's subtending edge
        all_edges.append(graft_recipient.edge)
    eligible_edges = [x for x in all_edges if filter_fn(x)]

    if not eligible_edges:
        raise Exception("could not place node {} in clade {}".format(graft, graft_recipient))
    focal_node = random.choice([x.head_node for x in eligible_edges])
    seed_node = focal_node.parent_node
    sisters = focal_node.sibling_nodes()

    # pick a child edge and detach its corresponding node
    #
    # DendroPy's Node.remove_child() messes with the edge lengths.
    # But, Node.clear_child_nodes() simply cuts that bit of the tree out.
    seed_node.clear_child_nodes()

    # set the correct edge length on the grafted node and make the grafted
    # node a child of the seed node
    graft.edge.length = seed_node.age - graft.age
    if graft.edge.length < 0:
        raise Exception("negative branch length")
    sisters.append(graft)
    seed_node.set_child_nodes(sisters)

    # make the focal node a child of the grafted node and set edge length
    focal_node.edge.length = graft.age - focal_node.age
    if focal_node.edge.length < 0:
        raise Exception("negative branch length")
    graft.add_child(focal_node)

    # return the (potentially new) crown of the clade
    if graft_recipient.parent_node == graft:
        return graft
    return graft_recipient


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
        [x.annotations.add_new("creation_method", "create_clade") for x in tree.preorder_node_iter()]
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
    [x.annotations.add_new("creation_method", "create_clade") for x in tree.preorder_node_iter()]
    assert is_binary(tree.seed_node.child_nodes()[0])
    tree.set_edge_lengths_from_node_ages(error_on_negative_edge_lengths=True)
    # Lock the child of the seed node so that things can still attach to the stem of this new clade
    lock_clade(tree.seed_node.child_nodes()[0])
    if list(get_short_branches(tree.seed_node)):
        logger.warning("{} short branches detected".format(len(list(get_short_branches(tree.seed_node)))))
    return tree


def lock_clade(node):
    pre = count_locked(node)
    for edge in edge_iter(node):
        edge.label = "locked"
    post = count_locked(node)
    if pre != post:
        logger.debug("locking clade: {} => {}".format(pre, post))


def count_locked(node):
    sum([x.label == "locked" for x in edge_iter(node)])


def is_fully_locked(node):
    return all([x.label == "locked" for x in edge_iter(node)])


def get_min_age(node):
    try:
        return min([x.head_node.age for x in edge_iter(node) if x.label != "locked"])
    except ValueError:
        return 0.0


def process_node(backbone_tree, backbone_bitmask, all_possible_tips, taxon_node, min_ccp, default_birth, default_death, yule=False):
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
        logger.debug("MRCA: skipping unlabeled rank with {} species".format(len(species)))
        return
    all_bitmask = backbone_tree.taxon_namespace.taxa_bitmask(labels=species)
    extant_bitmask = all_bitmask & backbone_bitmask
    if extant_bitmask is None or extant_bitmask == 0:
        logger.debug("MRCA: {} not present in backbone".format(taxon))
        mrca_rates[taxon] = (birth, death, 0.0, "from {} (unsampled)".format(parent))
        return
    mrca = backbone_tree.mrca(leafset_bitmask=extant_bitmask)
    if not species.issuperset(get_tip_labels(mrca)):
        logger.debug("MRCA: {} not monophyletic in backbone".format(taxon))
        mrca_rates[taxon] = (birth, death, 0.0, "from {} (not monophyletic)".format(parent))
        return
    extant = len(mrca.leaf_nodes())
    total = len(taxon_node.leaf_nodes())
    if extant > total:
        logger.warning("MRCA: {} has {} extant species but should have {} total species".format(taxon, extant, total))
        mrca_rates[taxon] = (birth, death, 0, "from {} (extant exceeds total)".format(parent))
        return
    ccp = crown_capture_probability(total, extant)
    if total == 1:
        logger.debug("MRCA: {} is a singleton".format(taxon))
        mrca_rates[taxon] = (birth, death, ccp, "from {} (singleton)".format(parent))
        return
    if ccp < min_ccp:
        logger.debug("MRCA: {} has crown capture probability {:.2f} < {:.2f} ({}/{} species)".format(taxon, ccp, min_ccp, extant, total))
        mrca_rates[taxon] = (birth, death, ccp, "from {} (crown capture probability)".format(parent))
        return
    sf = extant / total
    birth, death = get_birth_death_rates(mrca, sf, yule)
    logger.debug("MRCA: {} b={}, d={}, sf={}, ccp={}".format(taxon, birth, death, sf, ccp))
    mrca_rates[taxon] = (birth, death, ccp, "computed")


def run_precalcs(taxonomy_tree, backbone_tree, min_ccp=0.8, min_extant=3, yule=False):
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
    root_birth, root_death = get_birth_death_rates(root_mrca, len(root_mrca.leaf_nodes()) / len(all_possible_tips), yule)

    with click.progressbar(taxonomy_tree.preorder_internal_node_iter(exclude_seed_node=True), label="Rates", length=nnodes, show_pos=True, item_show_func=lambda x: x.label if x else None) as progress:
        for node in progress:
            # updates global mrca_rates as a side effect
            process_node(backbone_tree, backbone_bitmask, all_possible_tips, node, min_ccp, root_birth, root_death, yule)

    diff = time() - start_time
    if diff > 1:
        logger.debug("FastMRCA calculation time: {:.1f} seconds".format(diff))
    return mrca_rates


def update_tree_view(tree):
    # Stuff that DendroPy needs to keep a consistent view of the phylgoeny
    tree.calc_node_ages()
    tree.update_bipartitions()
    return get_tip_labels(tree)


def compute_node_depths(tree):
    res = dict()
    for leaf in tree.leaf_node_iter():
        cnt = 0
        for anc in leaf.ancestor_iter():
            if anc.label:
                cnt += 1
        res[leaf.taxon.label] = cnt
    return res


@click.command()
@click.option("--taxonomy", help="a taxonomy tree", type=click.File("r"), required=True)
@click.option("--backbone", help="the backbone tree to attach the taxonomy tree to", type=click.File("r"), required=True)
@click.option("--outgroups", help="comma separated list of outgroup taxa to ignore")
@click.option("--output", required=True, help="output base name to write out")
@click.option("--min-ccp", help="minimum probability to use to say that we've sampled the crown of a clade", default=0.8)
@click.option("--yule", help="assume a Yule pure-birth model (force extinction to be 0)", default=False, is_flag=True)
@click.option("-v", "--verbose", help="emit extra information (can be repeated)", count=True)
def main(taxonomy, backbone, outgroups, output, min_ccp, verbose, yule):
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

    logger.info("Reading taxonomy".format(taxonomy))
    taxonomy = dendropy.Tree.get_from_stream(taxonomy, schema="newick", rooting="default-rooted")
    tn = taxonomy.taxon_namespace
    tn.is_mutable = True
    if outgroups:
        outgroups = [x.replace("_", " ") for x in outgroups.split(",")]
        tn.new_taxa(outgroups)
    tn.is_mutable = False

    # Check for equal depth of all nodes
    node_depths = compute_node_depths(taxonomy)
    stats = collections.defaultdict(int)
    for v in node_depths.values():
        stats[v] += 1
    if len(stats) > 1:
        logger.warning("The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:")
        for k in sorted(stats.keys()):
            logger.warning("* {} tips have {} ranked ancestors".format(stats[k], k))
        logger.warning("If TACT-added tips are intruding into otherwise-monophyletic clades this should be corrected.")

    logger.info("Reading backbone".format(backbone))

    try:
        tree = dendropy.Tree.get_from_stream(backbone, schema="newick", rooting="default-rooted", taxon_namespace=tn)
    except dendropy.utility.error.ImmutableTaxonNamespaceError as e:
        logger.error("DendroPy error: {}".format(e.message))
        print("""
DendroPy error: {}

This usually indicates your backbone has species that are not present in your
taxonomy. Outgroups not in the taxonomy can be excluded with the argument:

    tact_add_taxa --outgroups outgroup_speciesA,outgroup_speciesB

For more details, run:

    tact_add_taxa --help
""".format(e.message))
        sys.exit(1)

    tree.encode_bipartitions()
    tree.calc_node_ages()

    tree_tips = get_tip_labels(tree)
    all_possible_tips = get_tip_labels(taxonomy)

    logger.info("Backbone needs to add {} tips".format(len(tree_tips.symmetric_difference(all_possible_tips))))

    full_clades = set()

    fastmrca.initialize(tree)

    with open(output + ".rates.csv", "w") as wfile:
        writer = csv.writer(wfile)
        writer.writerow(("taxon", "birth", "death", "ccp", "source"))
        for key, value in run_precalcs(taxonomy, tree, min_ccp, yule=yule).items():
            row = [key]
            row.extend(value)
            writer.writerow(row)

    initial_length = len(tree_tips)

    bar = click.progressbar(label="TACT",
                            length=len(all_possible_tips) - initial_length,
                            show_pos=True,
                            item_show_func=lambda x: x)

    def bar_update():
        bar.pos = len(tree_tips) - initial_length
        bar.current_item = taxon if taxon else ""
        bar.update(0)

    for taxon_node in taxonomy.postorder_internal_node_iter(exclude_seed_node=True):
        taxon = taxon_node.label
        if not taxon:
            continue
        species = get_tip_labels(taxon_node)
        extant_species = tree_tips.intersection(species)
        logger.info("**  {} ({}/{})  **".format(taxon, len(extant_species), len(species)))
        ccp = mrca_rates[taxon][2]

        clades_to_generate = full_clades.intersection([x.label for x in taxon_node.postorder_internal_node_iter(exclude_seed_node=True)])

        if not extant_species:
            # No species sampled, so create a clade from whole cloth
            logger.debug("    {}: no species sampled, will create later".format(taxon))
            full_clades.add(taxon)
            bar_update()
            continue

        # Check for monophyly for this node
        node = fastmrca.get(extant_species)
        if not node:
            logger.info("    {}: is not monophyletic".format(taxon))
            continue

        if extant_species == species:
            # Everything sampled and monophyletic, so skip this
            lock_clade(node)
            logger.debug("    {}: all species accounted for".format(taxon))
            continue

        if tree_tips.issuperset(species):
            # XXX: Does this check ever get triggered?
            lock_clade(node)
            logger.info("    {}: all species already present in tree".format(taxon))
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
                logger.info("    {}: skipping clade {} as all species already present in tree".format(taxon, clade))
                full_clades.remove(clade)
                continue
            logger.info("    {}: adding clade {} (n={})".format(taxon, clade, len(full_node.leaf_nodes())))
            # Generate all times needed to attach to the main clade
            times = get_new_branching_times(node, taxon_node, tree, tyoung=0, min_ccp=min_ccp, num_new_times=len(full_node_species))

            if is_fully_locked(node):
                logger.info("    {}: is fully locked, so attaching to stem".format(taxon))
                # Must attach to stem for this clade, so generate a time on the stem lineage
                times2 = get_new_branching_times(node, taxon_node, tree, min_ccp=min_ccp, told=node.parent_node.age, tyoung=node.age, num_new_times=1)
                # Drop the oldest time and add on our new time on the stem lineage
                times.sort()
                times.pop()
                times.append(times2.pop())
            else:
                # Even if the main clade isn't fully locked, it might have a constraint on a valid attachment point
                min_age = get_min_age(node)
                if min_age > 0 and max(times) < min_age:
                    logger.info("    {}: has a minimum age constraint {:.2f} but oldest generated time was {:.2f}".format(taxon, min_age, max(times)))
                    times2 = get_new_branching_times(node, taxon_node, tree, tyoung=min_age, min_ccp=min_ccp, num_new_times=1)
                    # Drop the oldest time and add on our new time on the stem lineage
                    times.sort()
                    times.pop()
                    times.append(times2.pop())

            # Generate a new tree
            new_tree = create_clade(tn, full_node_species, times)
            # Update our current MRCA node (because we might have attached to stem)
            node = graft_node(node, new_tree.seed_node, is_fully_locked(node) or ccp < min_ccp)
            # Stuff that DendroPy needs to keep a consistent view of the phylgoeny
            tree.calc_node_ages()
            tree.update_bipartitions()
            tree_tips = get_tip_labels(tree)
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
        logger.info("    {}: adding {} new species".format(taxon, len(species.difference(extant_species))))
        node = fastmrca.get(extant_species)
        times = get_new_branching_times(node, taxon_node, tree, tyoung=get_min_age(node), min_ccp=min_ccp)
        node = fill_new_taxa(tn, node, species.difference(tree_tips), times, ccp < min_ccp)
        # Update stuff
        tree_tips = update_tree_view(tree)
        # Since only monophyletic nodes get to here, lock this clade
        lock_clade(node)
        if not is_binary(node):
            # Shouldn't happen
            raise ValueError("Tree is not binary!")
        bar_update()

    assert(is_binary(tree.seed_node))
    # Reset terminal because we aren't using the context manager
    bar.render_finish()
    tree.ladderize()
    tree.write(path=output + ".newick.tre", schema="newick", suppress_rooting=True)
    tree.write(path=output + ".nexus.tre", schema="nexus")
    print()


if __name__ == '__main__':
    main()
