#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Try to assign tips to a pre-existing tree based on a taxonomy
# Jonathan Chang, May 13, 2016

from __future__ import division

import csv
import itertools
import collections
import sys
import random
from math import log, exp, ceil
from decimal import Decimal as D
from multiprocessing import Pool
import operator

import dendropy
import click

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)

from .lib import optim_bd, is_binary, get_short_branches, get_tip_labels, get_monophyletic_node, crown_capture_probability, edge_iter, get_new_times

def get_birth_death_rates(node, ages, sampling):
    logger.debug()

def get_new_branching_times(node, n_extant, n_total, told=None, tyoung=0, min_ccp=0.8):
    """
    Get `n_total` new branching times for a `node`.
    """
    if n_extant == n_total:
        raise Exception("get_new_branching_times args 2 and 3 cannot be equal")
    ccp = crown_capture_probability(n_total, n_extant)
    if n_extant == 1:
        # if we have a singleton then go up a node to get a better handle on
        # birth/death rates and origination times
        node = node.parent_node
        diff = n_total - n_extant
        n_extant = len(node.leaf_nodes())
        n_total = n_extant + diff
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    ages += [node.age]
    sampling = n_extant / n_total
    birth, death = get_birth_death_rates(node, ages, sampling)
    if node.annotations.get_value("birth"):
        birth = node.annotations.get_value("birth")
        death = node.annotations.get_value("death")
    else:
        birth, death = optim_bd(ages, sampling)
        node.annotations.add_new("birth", birth)
        node.annotations.add_new("death", death)
    if ccp < min_ccp and told is not None:
        told = node.parent_node.age
    return birth, death, ccp, get_new_times(ages, birth, death, n_total - n_extant, told, tyoung)


def get_new_branching_times(node, n_extant, n_total, told=None, tyoung=0, min_ccp=0.8):
    """
    Get `n_total` new branching times for a `node`.
    """
    if n_extant == n_total:
        raise Exception("get_new_branching_times args 2 and 3 cannot be equal")
    ccp = crown_capture_probability(n_total, n_extant)
    if n_extant == 1:
        # if we have a singleton then go up a node to get a better handle on
        # birth/death rates and origination times
        node = node.parent_node
        diff = n_total - n_extant
        n_extant = len(node.leaf_nodes())
        n_total = n_extant + diff
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    ages += [node.age]
    sampling = n_extant / n_total
    birth, death = get_birth_death_rates(node, ages, sampling)
    if node.annotations.get_value("birth"):
        birth = node.annotations.get_value("birth")
        death = node.annotations.get_value("death")
    else:
        birth, death = optim_bd(ages, sampling)
        node.annotations.add_new("birth", birth)
        node.annotations.add_new("death", death)
    if ccp < min_ccp and told is not None:
        told = node.parent_node.age
    return birth, death, ccp, get_new_times(ages, birth, death, n_total - n_extant, told, tyoung)

def fill_new_taxa(namespace, node, new_taxa, times, stem=False, excluded_nodes=None):
    if stem:
        node = node.parent_node

    for new_species, new_age in itertools.izip(new_taxa, times):
        new_node = dendropy.Node()
        new_node.annotations.add_new("creation_method", "fill_new_taxa")
        new_node.age = new_age
        new_leaf = new_node.new_child(taxon=namespace.require_taxon(new_species), edge_length=new_age)
        new_leaf.age = 0
        node = graft_node(node, new_node, stem)

    if list(get_short_branches(node)):
        logger.warn("{} short branches detected".format(len(list(get_short_branches(node)))))
        #import pdb; pdb.set_trace()

    node.locked = None

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
        node = random.sample(valid_nodes, 1).pop()
        child = node.new_child()
        child.age = age
    species = list(species)
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
        #import pdb; pdb.set_trace()
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

@click.command()
@click.option("--taxonomy", help="a taxonomy tree", type=click.File("rb"), required=True)
@click.option("--backbone", help="the backbone tree to attach the taxonomy tree to", type=click.File("rb"), required=True)
@click.option("--outgroups", help="comma separated list of outgroup taxa to ignore")
@click.option("--output", required=True, help="output base name to write out")
@click.option("--min-ccp", help="minimum probability to use to say that we've sampled the crown of a clade", default=0.8)
def main(taxonomy, backbone, outgroups, output, min_ccp):
    """
    Add tips onto a BACKBONE phylogeny using a TAXONOMY phylogeny.
    """
    logger.info("reading taxonomy")
    taxonomy = dendropy.Tree.get_from_stream(taxonomy, schema="newick")
    tn = taxonomy.taxon_namespace
    tn.is_mutable = True
    if outgroups:
        outgroups = [x.replace("_", " ") for x in outgroups.split(",")]
        tn.new_taxa(outgroups)
    tn.is_mutable = False

    logger.info("reading tree")
    tree = dendropy.Tree.get_from_stream(backbone, schema="newick", rooting="force-rooted", taxon_namespace=tn)
    tree.encode_bipartitions()
    tree.calc_node_ages()

    tree_tips = get_tip_labels(tree)
    all_possible_tips = get_tip_labels(taxonomy)

    logger.info("{} tips to add".format(len(tree_tips.symmetric_difference(all_possible_tips))))

    full_clades = set()

    for taxon_node in taxonomy.postorder_internal_node_iter(exclude_seed_node=True):
        taxon = taxon_node.label
        if not taxon:
            continue
        spaces = taxon_node.level() * "  "
        species = set([x.taxon.label for x in taxon_node.leaf_iter()])
        extant_species = tree_tips.intersection(species)
        logger.info("{}{} ({}/{})... ({} remain)".format(spaces, taxon, len(extant_species), len(species), len(tree_tips.symmetric_difference(all_possible_tips))))

        clades_to_generate = full_clades.intersection([x.label for x in taxon_node.postorder_internal_node_iter(exclude_seed_node=True)])
        to_remove = set([])

        if extant_species:
            if extant_species == species:
                #logger.debug("{}  => all species accounted for".format(spaces))
                continue
            if tree_tips.issuperset(species):
                logger.info("{}  => all species already present in tree".format(spaces))
                continue

            node = get_monophyletic_node(tree, extant_species)
            if not node:
                logger.info(spaces + "  => not monophyletic")
                continue

            clade_sizes = [(clade, len(taxonomy.find_node_with_label(clade).leaf_nodes())) for clade in clades_to_generate]

            # sorting clades by size should add genera before families... better way would be to sort by rank
            for clade, clade_size in sorted(clade_sizes, key=operator.itemgetter(1)):
                full_node = taxonomy.find_node_with_label(clade)
                full_node_species = [x.taxon.label for x in full_node.leaf_iter()]
                if tree_tips.issuperset(full_node_species):
                    logger.info("{}  => skipping {} as all species already present in tree".format(spaces, clade))
                    full_clades.remove(clade)
                    continue
                birth, death, ccp, times = get_new_branching_times(node, len(species), len(species) + len(full_node_species), tyoung=get_min_age(node), min_ccp=min_ccp)
                logger.info("{}  => adding {} (n={})".format(spaces, clade, clade_size))

                if is_fully_locked(node):
                    # must attach to stem for this clade, so generate a time on the stem lineage
                    _, _, _, times2 = get_new_branching_times(node, len(species), len(species) + 1, told=node.parent_node.age, tyoung=node.age, min_ccp=min_ccp)
                    times.sort()
                    times.pop()
                    times.append(times2.pop())

                # generate a new tree
                new_tree = create_clade(tn, full_node_species, times)
                node = graft_node(node, new_tree.seed_node, is_fully_locked(node))
                tree.calc_node_ages()
                tree.update_bipartitions()
                tree_tips = get_tip_labels(tree)
                extant_species = tree_tips.intersection(species)
                full_clades.remove(clade)
                assert(is_binary(node))

            # check to see if we need to continue adding species
            if extant_species == species:
                # lock clade since it is monophyletic and filled
                lock_clade(node)
                continue
            if len(extant_species) == len(species):
                raise Exception("enough species are present but mismatched?")

            logger.info("{}  => adding {} new species".format(spaces, len(species.difference(extant_species))))
            node = get_monophyletic_node(tree, extant_species)
            birth, death, ccp, times = get_new_branching_times(node, len(extant_species), len(species), tyoung=get_min_age(node), min_ccp=min_ccp)
            fill_new_taxa(tn, node, species.difference(tree_tips), times, ccp < min_ccp)
            tree.update_bipartitions()
            tree.calc_node_ages()
            tree_tips = get_tip_labels(tree)
            # since only monophyletic nodes get to here, lock this clade
            lock_clade(node)
            assert(is_binary(node))
        else:
            # create clade from whole cloth
            full_clades.add(taxon)

    assert(is_binary(tree.seed_node))
    tree.ladderize()
    for leaf in tree.leaf_node_iter():
        if leaf.edge.length <= 0.001:
            logger.info("warning: taxon {} has extremely short branch ({})".format(leaf.taxon.label, leaf.edge.length))
    tree.write(path=output + ".newick.tre", schema="newick")
    tree.write(path=output + ".nexus.tre", schema="nexus")

if __name__ == '__main__':
    main()

"""
tmp = tree.extract_tree_with_taxa([x.taxon for x in node.leaf_iter()])
tmp.write_to_path("tmp.tre", schema="newick")
"""
