#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Try to assign tips to a pre-existing tree based on a taxonomy
# Jonathan Chang, May 13, 2016

# Pragmata
from __future__ import division

# Python standard library
import Queue
import functools
import itertools
import logging
import operator
import random
import sys
from time import time
import math
import collections

logger = logging.getLogger(__name__)

# Third party
import dendropy
import click

# Internal
from .lib import optim_bd, is_binary, get_short_branches, get_tip_labels, get_monophyletic_node, crown_capture_probability, edge_iter, get_new_times
from . import fastmrca

global invalid_map
invalid_map = {}
def search_ancestors_for_valid_backbone_node(taxonomy_node, backbone_tips, ccp):
    global invalid_map
    seen = []
    target_node = None
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
    seen.pop() # ignore last node
    for x in seen:
        invalid_map[x] = taxonomy_target
    return (taxonomy_target, backbone_target)

def get_birth_death_rates(node, sampfrac):
    if len(node.child_nodes()) == 0: # is a tip
        return [math.log(1. / sampfrac) / node.parent_node.age, 0] # M-S estimate of rates
    return optim_bd(get_ages(node), sampfrac)

def get_ages(node):
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    ages += [node.age]
    return ages

def get_new_branching_times(backbone_node, taxonomy_node, backbone_tree, told=None, tyoung=0, min_ccp=0.8, num_new_times=None):
    """
    Get `n_total` new branching times for a `node`.
    """
    original_backbone_node = backbone_node
    original_taxonomy_node = taxonomy_node
    n_extant = len(original_backbone_node.leaf_nodes())
    n_total = len(original_taxonomy_node.leaf_nodes())
    if num_new_times is None:
        num_new_times = n_total - n_extant
    new_ccp = ccp = crown_capture_probability(n_total, n_extant)
    if new_ccp < min_ccp:
        full_tax = get_tip_labels(taxonomy_node)
        extant_tax = full_tax.intersection(get_tip_labels(backbone_node))
        mrca_node = fastmrca.get(extant_tax)
        if mrca_node is not None:
            logger.info("    {}: is monophyletic but has poor sampling; using stem (ccp {:.2f} < min_ccp {})".format(taxonomy_node.label, new_ccp, min_ccp))
            backbone_node = backbone_node.parent_node
        else:
            while new_ccp < min_ccp:
                logger.info("    {}: is not monophyletic and has poor sampling, checking ancestors (ccp {:.2f} < min_ccp {})".format(taxonomy_node.label, new_ccp, min_ccp))
                taxonomy_node, backbone_node = search_ancestors_for_valid_backbone_node(taxonomy_node, get_tip_labels(backbone_tree), ccp=min_ccp)
                n_extant = len(backbone_node.leaf_nodes())
                n_total = len(taxonomy_node.leaf_nodes())
                new_ccp = crown_capture_probability(n_total, n_extant)
    sampling = n_extant / n_total
    if backbone_node.annotations.get_value("birth"):
        #logger.debug("    {}: Cache HIT on birth/death rates".format(taxonomy_node.label))
        birth = backbone_node.annotations.get_value("birth")
        death = backbone_node.annotations.get_value("death")
    else:
        logger.debug("    {}: Cache MISS on birth/death rates".format(taxonomy_node.label))
        birth, death = get_birth_death_rates(backbone_node, sampling)
        backbone_node.annotations.add_new("birth", birth)
        backbone_node.annotations.add_new("death", death)
    if ccp < min_ccp and told is not None:
        told = original_backbone_node.parent_node.age
    if len(original_backbone_node.leaf_nodes()) == 1 and told is None:
        # attach to stem in the case of a singleton
        told = original_backbone_node.parent_node.age
    times = get_new_times(get_ages(original_backbone_node), birth, death, num_new_times, told, tyoung)
    return birth, death, ccp, times

def fill_new_taxa(namespace, node, new_taxa, times, stem=False, excluded_nodes=None):
    # lol, graft_node already accounts for this so don't do it here!!
    #if stem:
        #node = node.parent_node

    for new_species, new_age in itertools.izip(new_taxa, times):
        new_node = dendropy.Node()
        new_node.annotations.add_new("creation_method", "fill_new_taxa")
        new_node.age = new_age
        new_leaf = new_node.new_child(taxon=namespace.require_taxon(new_species), edge_length=new_age)
        new_leaf.annotations.add_new("creation_method", "fill_new_taxa")
        new_leaf.age = 0
        node = graft_node(node, new_node, stem)

    if list(get_short_branches(node)):
        logger.warn("{} short branches detected".format(len(list(get_short_branches(node)))))

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
    lock_clade(tree.seed_node)
    if list(get_short_branches(tree.seed_node)):
        logger.warn("{} short branches detected".format(len(list(get_short_branches(tree.seed_node)))))
    return tree

def lock_clade(node):
    for edge in edge_iter(node):
        edge.label = "locked"

def is_fully_locked(node):
    return all([x.label == "locked" for x in edge_iter(node)])

def get_min_age(node):
    try:
        return min([x.head_node.age for x in edge_iter(node) if x.label is not "locked"])
    except ValueError:
        return 0.0

def process_node(backbone_tree, backbone_bitmask, all_possible_tips, taxon_node):
    taxon = taxon_node.label
    if not taxon:
        # ignore unlabeled ranks
        return None
    species = get_tip_labels(taxon_node)
    all_bitmask = backbone_tree.taxon_namespace.taxa_bitmask(labels=species)
    extant_bitmask = all_bitmask & backbone_bitmask
    if extant_bitmask is None or extant_bitmask == 0:
        # clade doesn't exist?
        return (taxon_node, None, None, None)
    mrca = backbone_tree.mrca(leafset_bitmask=extant_bitmask)
    if mrca:
        birth, death = get_birth_death_rates(mrca, len(mrca.leaf_nodes()) / len(taxon_node.leaf_nodes()))
        return (taxon_node, all_bitmask, birth, death)
    else:
        return (taxon_node, None, None, None)

def run_precalcs(taxonomy_tree, backbone_tree, min_ccp=0.8, min_extant=3):
    tree_tips = get_tip_labels(backbone_tree)
    backbone_bitmask = fastmrca.bitmask(tree_tips)
    all_possible_tips = get_tip_labels(taxonomy_tree)

    nnodes = len(taxonomy_tree.internal_nodes(exclude_seed_node=True))

    start_time = time()

    def annotate_result_node(result):
        if result is None:
            return
        taxon_node, taxon_bitmask, birth, death = result
        if birth is None or death is None:
            return
        backbone_node = backbone_tree.mrca(leafset_bitmask=taxon_bitmask & backbone_bitmask)
        if backbone_node:
            backbone_node.annotations.add_new("birth", birth)
            backbone_node.annotations.add_new("death", death)

    if fastmrca.cores == 1 or nnodes < 500:
        logger.debug("Precomputing rates serially since cores=1 ({}) or nnodes < 500 ({})".format(fastmrca.cores, nnodes))
        with click.progressbar(taxonomy_tree.preorder_internal_node_iter(exclude_seed_node=True), label="Calculating rates", length=nnodes, show_pos=True, item_show_func=lambda x: x.label if x else None) as progress:
            for node in progress:
                annotate_result_node(process_node(backbone_tree, backbone_bitmask, all_possible_tips, node))
    else:
        # While it would be incredibly easy to just run pool.unordered_imap
        # on everything, in practice this doesn't work because of the huge
        # variation in the time it takes to compute the birth/death rates
        # based on the number of tips in each subclade. So bucket each node
        # that we have to calculate the rate for based on the number of tips
        # that descend from that node, such that each bucket has roughly
        # the same number of tips (not nodes!)
        buckets = []
        tips_per_node = [(len(x.leaf_nodes()), x) for x in taxonomy_tree.preorder_internal_node_iter(exclude_seed_node=True)]
        queue = Queue.PriorityQueue()
        for x in range(max(int(fastmrca.cores/4), 2)):
            buckets.append([])
            queue.put((0, x))
        sums = [0] * fastmrca.cores

        for ntips, node in sorted(tips_per_node, key=operator.itemgetter(1), reverse=True):
            _, i = queue.get()
            buckets[i].append(node)
            sums[i] += ntips
            queue.put((sums[i], i))

        buckets.sort(key=lambda x: len(x))
        logger.debug("Precomputing rates in parallel, worker assignments ({} cores): {}".format(len(buckets), [len(x) for x in buckets]))

        progress = click.progressbar(label="Calculating rates", length=nnodes, show_pos=True)

        # Submit to the pool and keep track of promises...
        promises = []
        fn = functools.partial(process_node, backbone_tree, backbone_bitmask, all_possible_tips)
        for acc_nodes in buckets:
            promises.append(fastmrca.pool.map_async(fn, acc_nodes, len(acc_nodes)))

        # Ideally we would asynchronously resolve promises but this
        # doesn't work for some reason, so just resolve them in order
        # despite the worse UX
        for promise in promises:
            results = promise.get()
            for result in results:
                progress.update(1)
                annotate_result_node(result)

    diff = time() - start_time
    if diff > 1:
        logger.debug("FastMRCA calculation time: {:.1f} seconds".format(diff))

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
@click.option("--taxonomy", help="a taxonomy tree", type=click.File("rb"), required=True)
@click.option("--backbone", help="the backbone tree to attach the taxonomy tree to", type=click.File("rb"), required=True)
@click.option("--outgroups", help="comma separated list of outgroup taxa to ignore")
@click.option("--output", required=True, help="output base name to write out")
@click.option("--min-ccp", help="minimum probability to use to say that we've sampled the crown of a clade", default=0.8)
@click.option("--cores", help="maximum number of cores to use for parallel operations", type=int)
@click.option("-v", "--verbose", help="emit extra information (can be repeated)", count=True)
@click.option("--log", "log_file", help="if verbose output is enabled, send it to this file instead of standard output")
def main(taxonomy, backbone, outgroups, output, min_ccp, cores, verbose, log_file):
    """
    Add tips onto a BACKBONE phylogeny using a TAXONOMY phylogeny.
    """

    if verbose >= 2:
        logger.setLevel(logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)
    if log_file:
        logger.addHandler(logging.FileHandler(log_file))
    else:
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
    for v in node_depths.itervalues():
        stats[v] += 1
    if len(stats) > 1:
        logger.warn("The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:")
        for k in sorted(stats.keys()):
            logger.warn("* {} tips have {} ranked ancestors".format(stats[k], k))
        logger.warn("If TACT-added tips are intruding into otherwise-monophyletic clades this should be corrected.")

    logger.info("Reading backbone".format(backbone))

    try:
        tree = dendropy.Tree.get_from_stream(backbone, schema="newick", rooting="default-rooted", taxon_namespace=tn)
    except dendropy.utility.error.ImmutableTaxonNamespaceError as e:
        logger.error("DendroPy error: {}".format(e.message))
        print """
DendroPy error: {}

This usually indicates your backbone has species that are not present in your
taxonomy. Outgroups not in the taxonomy can be excluded with the argument:

    tact_add_taxa --outgroups outgroup_speciesA,outgroup_speciesB

For more details, run:

    tact_add_taxa --help
""".format(e.message)
        sys.exit(1)

    tree.encode_bipartitions()
    tree.calc_node_ages()

    tree_tips = get_tip_labels(tree)
    all_possible_tips = get_tip_labels(taxonomy)

    logger.info("Backbone needs to add {} tips".format(len(tree_tips.symmetric_difference(all_possible_tips))))

    full_clades = set()

    fastmrca.initialize(tree)
    logger.debug("FastMRCA autotuned parameters: single-thread cutoff is {}".format(fastmrca.maxtax))

    run_precalcs(taxonomy, tree, min_ccp)

    initial_length = len(tree_tips)

    bar = click.progressbar(label="Adding taxa",
            length=len(all_possible_tips) - initial_length,
            show_pos=True,
            item_show_func=lambda x: x)

    def bar_update():
        bar.pos = len(tree_tips) - initial_length;
        bar.current_item = taxon if taxon else ""
        bar.update(0)

    for taxon_node in taxonomy.postorder_internal_node_iter(exclude_seed_node=True):
        taxon = taxon_node.label
        if not taxon:
            continue
        species = get_tip_labels(taxon_node)
        extant_species = tree_tips.intersection(species)
        logger.info("**  {} ({}/{})  **".format(taxon, len(extant_species), len(species)))

        clades_to_generate = full_clades.intersection([x.label for x in taxon_node.postorder_internal_node_iter(exclude_seed_node=True)])
        to_remove = set([])

        if not extant_species:
            # No species sampled, so create a clade from whole cloth
            logger.debug("    {}: no species sampled, will create later".format(taxon))
            full_clades.add(taxon)
            bar_update()
            continue

        if extant_species == species:
            # Everything sampled, so skip this
            logger.debug("    {}: all species accounted for".format(taxon))
            continue

        if tree_tips.issuperset(species):
            # XXX: Does this check ever get triggered?
            logger.info("    {}: all species already present in tree".format(taxon))
            continue

        node = fastmrca.get(extant_species)
        if not node:
            logger.info("    {}: is not monophyletic".format(taxon))
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
            birth, death, ccp, times = get_new_branching_times(node, taxon_node, tree, tyoung=get_min_age(node), min_ccp=min_ccp, num_new_times=len(full_node_species))
            logger.info("    {}: adding clade {} (n={})".format(taxon, clade, len(full_node.leaf_nodes())))

            if is_fully_locked(node):
                logger.info("    {}: clade {} is fully locked, so attaching to stem".format(taxon, taxon))
                # Must attach to stem for this clade, so generate a time on the stem lineage
                _, _, _, times2 = get_new_branching_times(node, taxon_node, tree, min_ccp=min_ccp, told=node.parent_node.age, tyoung=node.age, num_new_times=1)
                # Drop the oldest time and add on our new time on the stem lineage
                times.sort()
                times.pop()
                times.append(times2.pop())

            # Generate a new tree
            new_tree = create_clade(tn, full_node_species, times)
            new_node = graft_node(node, new_tree.seed_node, is_fully_locked(node))
            # Stuff that DendroPy needs to keep a consistent view of the phylgoeny
            tree.calc_node_ages()
            tree.update_bipartitions()
            tree_tips = get_tip_labels(tree)
            # Update our view of what's in the tree
            extant_species = tree_tips.intersection(species)
            # We've added this clade so pop it off our stack
            full_clades.remove(clade)
            if not is_binary(new_node):
                raise ValueError("New grafted node is not binary!")
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
        birth, death, ccp, times = get_new_branching_times(node, taxon_node, tree, tyoung=get_min_age(node), min_ccp=min_ccp)
        fill_new_taxa(tn, node, species.difference(tree_tips), times, ccp < min_ccp)
        # Update stuff
        tree_tips = update_tree_view(tree)
        # Reacquire the potentially new MRCA of this clade, with everything added
        node = fastmrca.get(species)
        # Since only monophyletic nodes get to here, lock this clade
        lock_clade(node)
        if not is_binary(node):
            # Shouldn't happen
            raise ValueError("Tree is not binary!")
        bar_update()

    fastmrca.cleanup()
    assert(is_binary(tree.seed_node))
    tree.ladderize()
    tree.write(path=output + ".newick.tre", schema="newick")
    tree.write(path=output + ".nexus.tre", schema="nexus")
    print

if __name__ == '__main__':
    main()

"""
tmp = tree.extract_tree_with_taxa([x.taxon for x in node.leaf_iter()])
tmp.write_to_path("tmp.tre", schema="newick")
"""
