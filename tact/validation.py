# -*- coding: utf-8 -*-

"""Various validation functions for inputs."""

import collections

import click
import dendropy

from .tree_util import compute_node_depths


def validate_outgroups(ctx, param, value):
    if value is None:
        return
    try:
        value = value.split(",")
    except AttributeError:
        # Tuples and lists shouldn't have the .split method
        pass
    return [x.replace("_", " ") for x in value]


def validate_newick(ctx, param, value):
    return dendropy.Tree.get_from_stream(value, schema="newick", rooting="default-rooted")


def validate_tree_node_depths(ctx, param, value):
    node_depths = compute_node_depths(value)
    stats = collections.defaultdict(int)
    for v in node_depths.values():
        stats[v] += 1
    if len(stats) > 1:
        msg = "The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:\n"
        for k in sorted(stats.keys()):
            msg += f"* {stats[k]} tips have {k} ranked ancestors\n"
        raise click.BadParameter(msg)
    return value


def validate_taxonomy_tree(ctx, param, value):
    value = validate_newick(ctx, param, value)
    return validate_tree_node_depths(ctx, param, value)
