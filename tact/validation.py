"""Various validation functions for `click` classes and parameters."""

import collections
from collections.abc import Sequence
from typing import IO, Any

import click
import dendropy
from dendropy.utility.error import ImmutableTaxonNamespaceError

from .cli_taxonomy import echo_changes, read_taxonomy
from .tree_util import compute_node_depths, is_binary, is_ultrametric, update_tree_view


def validate_outgroups(
    ctx: click.Context | None, param: click.Parameter | None, value: str | Sequence[str] | None
) -> list[str] | None:
    """Validates an `outgroups` parameter, by splitting on commas and transforming underscores to spaces."""
    if value is None:
        return None
    if isinstance(value, str):
        value = value.split(",")
    return [x.replace("_", " ") for x in value]


def validate_newick(
    ctx: click.Context | None, param: click.Parameter | None, value: IO[str], **kwargs: object
) -> dendropy.Tree:
    """Validates a Newick tree, using appropriate defaults."""
    return dendropy.Tree.get_from_stream(value, schema="newick", rooting="default-rooted", **kwargs)


def validate_tree_node_depths(
    ctx: click.Context | None, param: click.Parameter | None, value: dendropy.Tree
) -> dendropy.Tree:
    """Validates a DendroPy tree, ensuring that the node depth is equal for all tips."""
    node_depths = compute_node_depths(value)
    stats: collections.defaultdict[int, int] = collections.defaultdict(int)
    for v in node_depths.values():
        stats[v] += 1
    if len(stats) > 1:
        msg = "The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:\n"
        for k in sorted(stats.keys()):
            msg += f"* {stats[k]} tips have {k} ranked ancestors\n"
        raise click.BadParameter(msg)
    return value


def validate_taxonomy_tree(ctx: click.Context | None, param: click.Parameter | None, value: str) -> dendropy.Tree:
    """Validates a taxonomy, given either as a Newick tree or as a CSV file that is built into a tree in memory.

    The format is sniffed from the first character of the file rather than taken from its extension, because
    taxonomy trees go by many extensions (.tre, .tree, .nwk, .txt) and a CSV need not end in .csv. A Newick
    tree starts with `(`, or `[` for a comment such as `[&R]`, which a CSV header row does not. Anything else
    is read as a CSV, so that a file that is neither gets the list of problems from `read_taxonomy` rather than
    a DendroPy parse error.

    If the command has a `fill_empty` parameter, it must be eager so that it is set before this callback runs.
    """
    fill_empty = bool(ctx and ctx.params.get("fill_empty"))
    # utf-8-sig drops the byte order mark that spreadsheet programs often add
    with open(value, encoding="utf-8-sig") as fh:
        text = fh.read()
    start = text.lstrip()
    if start[:6].upper() == "#NEXUS":
        raise click.BadParameter("this is a Nexus file; give a Newick tree or a CSV file instead")
    if start[:1] in ("(", "["):
        if fill_empty:
            raise click.BadParameter("--fill-empty only applies to a taxonomy CSV, not a tree")
        tree = dendropy.Tree.get(data=text, schema="newick", rooting="default-rooted")
    else:
        taxonomy = read_taxonomy(value, fill_empty)
        echo_changes(taxonomy)
        tree = dendropy.Tree.get(data=taxonomy.to_newick(), schema="newick", rooting="default-rooted")
    return validate_tree_node_depths(ctx, param, tree)


class BackboneCommand(click.Command):
    """Helper class to validate a Click Command that contains a backbone tree.

    At a minimum, the Command must contain a `backbone` parameter, which is validated by `validate_newick`
    and checked to ensure it is a binary tree.

    If the command also contains a `taxonomy` parameter, representing a taxonomic phylogeny,
    this is also validated to ensure that the DendroPy TaxonNamespace is non-strict superset
    of the taxa contained in `backbone`. An optional `outgroups` parameter may add
    other taxa not in the `taxonomy`.

    If the command also contains an `ultrametricity_precision` parameter, the
    ultrametricity of the `backbone` is also checked.
    """

    def validate_backbone_variables(self, ctx: click.Context, params: dict[str, Any]) -> dict[str, Any]:
        """Validates variables related to the backbone and taxonomy files."""
        if "taxonomy" in params:
            tn = params["taxonomy"].taxon_namespace
            tn.is_mutable = True
            if params.get("outgroups"):
                tn.new_taxa(params["outgroups"])
            tn.is_mutable = False
            try:
                backbone = validate_newick(ctx, None, params["backbone"], taxon_namespace=tn)
            except ImmutableTaxonNamespaceError as e:
                msg = f"""
                DendroPy error: {e}

                This usually indicates your backbone has species that are not present in your
                taxonomy. Outgroups not in the taxonomy can be excluded with the --outgroups argument.
                """
                raise click.BadParameter(msg) from None
        else:
            backbone = validate_newick(ctx, None, params["backbone"])

        if not is_binary(backbone.seed_node):
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
                raise click.BadParameter(msg) from None

        params["backbone"] = backbone
        return params

    def make_context(
        self, info_name: str | None, args: list[str], parent: click.Context | None = None, **extra: object
    ) -> click.Context:
        """Set up the proper Click context for a command handler."""
        ctx = super().make_context(info_name, args, parent, **extra)
        ctx.params = self.validate_backbone_variables(ctx, ctx.params)
        return ctx
