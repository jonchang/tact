# -*- coding: utf-8 -*-

"""Various validation functions for `click` classes and parameters."""

import collections

import click
import dendropy

from .tree_util import compute_node_depths
from .tree_util import is_binary
from .tree_util import is_ultrametric
from .tree_util import update_tree_view


def validate_outgroups(ctx, param, value):
    """Validates an `outgroups` parameter, by splitting on commas and transforming underscores to spaces."""
    if value is None:
        return
    try:
        value = value.split(",")
    except AttributeError:
        # Tuples and lists shouldn't have the .split method
        pass
    return [x.replace("_", " ") for x in value]


def validate_newick(ctx, param, value, **kwargs):
    """Validates a Newick tree, using appropriate defaults."""
    return dendropy.Tree.get_from_stream(value, schema="newick", rooting="default-rooted", **kwargs)


def validate_tree_node_depths(ctx, param, value):
    """Validates a DendroPy tree, ensuring that the node depth is equal for all tips."""
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
    """Validates a taxonomy tree."""
    value = validate_newick(ctx, param, value)
    return validate_tree_node_depths(ctx, param, value)


class BackboneCommand(click.Command):
    """
    Helper class to validate a Click Command that contains a backbone tree.

    At a minimum, the Command must contain a `backbone` parameter, which is validated by `validate_newick`
    and checked to ensure it is a binary tree.

    If the command also contains a `taxonomy` parameter, representing a taxonomic phylogeny,
    this is also validated to ensure that the DendroPy TaxonNamespace is non-strict superset
    of the taxa contained in `backbone`. An optional `outgroups` parameter may add
    other taxa not in the `taxonomy`.

    If the command also contains an `ultrametricity_precision` parameter, the
    ultrametricity of the `backbone` is also checked.
    """

    def validate_backbone_variables(self, ctx, params):
        if "taxonomy" in params:
            tn = params["taxonomy"].taxon_namespace
            tn.is_mutable = True
            if "outgroups" in params and params["outgroups"]:
                tn.new_taxa(params["outgroups"])
            tn.is_mutable = False
            try:
                backbone = validate_newick(ctx, params, params["backbone"], taxon_namespace=tn)
            except dendropy.utility.error.ImmutableTaxonNamespaceError as e:
                msg = f"""
                DendroPy error: {e}

                This usually indicates your backbone has species that are not present in your
                taxonomy. Outgroups not in the taxonomy can be excluded with the --outgroups argument.
                """
                raise click.BadParameter(msg)
        else:
            backbone = validate_newick(ctx, params, params["backbone"])

        if not is_binary(backbone):
            raise click.BadParameter("Backbone tree is not binary!")
        update_tree_view(backbone)

        if "ultrametricity_precision" in params:
            ultra, res = is_ultrametric(backbone, params["ultrametricity_precision"])
            if not ultra:
                msg = f"""
                Tree is not ultrametric!
                {res[0][0]} has a root distance of {res[0][1]}, but {res[1][0]} has {res[1][1]}

                Increase `--ultrametricity-precision` or use phytools::force.ultrametric in R
                """
                raise click.BadParameter(msg)

        params["backbone"] = backbone
        return params

    def make_context(self, *args, **kwargs):
        ctx = super(BackboneCommand, self).make_context(*args, **kwargs)
        ctx.params = self.validate_backbone_variables(ctx, ctx.params)
        return ctx
