"""Various validation functions for `click` classes and parameters."""

import collections
import math
import os
from collections.abc import Sequence
from typing import IO, Any

import click
import dendropy
from dendropy.dataio.newickreader import NewickReader
from dendropy.utility.error import DataParseError

from .cli_taxonomy import MAX_PROBLEMS, ROOT, echo_changes, read_taxonomy
from .names import adopt_names, clean_name, clean_tree_names, name_key
from .report import fmt_species_list
from .tree_util import compute_node_depths, is_ultrametric, update_tree_view

NEXUS_TREE = "this is a Nexus file! Use a Newick tree instead"


def is_nexus(text: str) -> bool:
    """Whether `text` is the start of a Nexus file."""
    return text.lstrip()[:6].upper() == "#NEXUS"


def fmt_spellings(spellings: Sequence[str]) -> str:
    """Join the spellings of a name, as `A, B and C`."""
    return ", ".join(spellings[:-1]) + " and " + spellings[-1]


def fmt_bullets(items: Sequence[str]) -> str:
    """List `items` one per line, up to `MAX_PROBLEMS` of them, saying how many more there are."""
    lines = [f"* {item}" for item in items[:MAX_PROBLEMS]]
    if len(items) > MAX_PROBLEMS:
        lines.append(f"* ... and {len(items) - MAX_PROBLEMS} more")
    return "\n".join(lines)


def validate_output(ctx: click.Context | None, param: click.Parameter | None, value: str) -> str:
    """Validates an output base name, checking that its directory exists so that writing to it later will work.

    Make the parameter eager if other parameters' callbacks are slow, so that this runs before them.
    """
    directory = os.path.dirname(value) or "."
    if not os.path.isdir(directory):
        raise click.BadParameter(f"directory '{directory}' does not exist")
    if not os.access(directory, os.W_OK):
        raise click.BadParameter(f"directory '{directory}' is not writable")
    return value


def validate_newick_path(ctx: click.Context | None, param: click.Parameter | None, value: str) -> str:
    """Validates that a path looks like a Newick tree, from its first characters, without reading the entire tree."""
    with open(value, encoding="utf-8-sig", errors="replace") as fh:
        start = fh.read(1024)
    if is_nexus(start):
        raise click.BadParameter(NEXUS_TREE)
    if start.lstrip()[:1] not in ("(", "["):
        msg = "this is not a Newick tree, which starts with '('"
        if param is not None and param.name == "taxonomy":
            msg += "; make one from a taxonomy CSV with tact_build_taxonomic_tree"
        raise click.BadParameter(msg)
    return value


def validate_outgroups(
    ctx: click.Context | None, param: click.Parameter | None, value: str | Sequence[str] | None
) -> list[str] | None:
    """Validates an `outgroups` parameter, by splitting on commas and cleaning each name."""
    if value is None:
        return None
    if isinstance(value, str):
        value = value.split(",")
    return [clean_name(x) for x in value]


def validate_newick(ctx: click.Context | None, param: click.Parameter | None, value: IO[str]) -> dendropy.Tree:
    """Validates a Newick tree, using appropriate defaults."""
    return read_newick(ctx, param, value.read())


def read_newick(ctx: click.Context | None, param: click.Parameter | None, text: str) -> dendropy.Tree:
    """Read a Newick tree into a namespace of its own, with its names cleaned.

    Names that differ only in spaces, underscores or case are an error, rather than being silently merged.
    """
    if is_nexus(text):
        raise click.BadParameter(NEXUS_TREE, ctx, param)
    if not text.strip():
        raise click.BadParameter("the file is empty", ctx, param)
    if not text.rstrip().endswith(";"):
        # DendroPy suggests one of its own options for this
        raise click.BadParameter("this is not a Newick tree: it does not end with a semicolon", ctx, param)
    try:
        # Underscores and case are kept so that the names are compared and cleaned in one place
        tree = dendropy.Tree.get(
            data=text,
            schema="newick",
            rooting="default-rooted",
            preserve_underscores=True,
            case_sensitive_taxon_labels=True,
            taxon_namespace=dendropy.TaxonNamespace(is_case_sensitive=True),
        )
    except NewickReader.NewickReaderDuplicateTaxonError as e:
        # DendroPy puts the label at the end of a long explanation
        label = clean_name((e.message or "").rsplit(": ", 1)[-1])
        raise click.BadParameter(f"the tip {label} appears more than once", ctx, param) from None
    except DataParseError as e:
        msg = f"this is not a Newick tree: {e.message}"
        if e.line_num is not None:
            msg += f" (line {e.line_num}, column {e.col_num})"
        raise click.BadParameter(msg, ctx, param) from None
    tip_clashes, taxon_clashes = clean_tree_names(tree, keep=(ROOT,))
    problems = [f"tips {fmt_spellings(group)}" for group in tip_clashes]
    problems += [f"taxa {fmt_spellings(group)}" for group in taxon_clashes]
    if problems:
        msg = (
            f"{len(problems)} name{'s are' if len(problems) != 1 else ' is'} spelled more than one way, "
            f"differing only in spaces, underscores or case:\n{fmt_bullets(problems)}"
        )
        raise click.BadParameter(msg, ctx, param)
    # Now that no two tips differ only in case, looking them up by label can ignore it
    tree.taxon_namespace.is_case_sensitive = False
    return tree


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
    if is_nexus(text):
        raise click.BadParameter("this is a Nexus file; give a Newick tree or a CSV file instead")
    if text.lstrip()[:1] in ("(", "["):
        if fill_empty:
            raise click.BadParameter("--fill-empty only applies to a taxonomy CSV, not a tree")
        tree = read_newick(ctx, param, text)
    else:
        taxonomy = read_taxonomy(value, fill_empty)
        echo_changes(taxonomy)
        tree = read_newick(ctx, param, taxonomy.to_newick())
    return validate_tree_node_depths(ctx, param, tree)


def first_tip(node: dendropy.Node) -> str:
    """The label of the first tip below `node`."""
    return next(node.leaf_iter()).taxon.label


def check_binary(ctx: click.Context | None, param: click.Parameter | None, tree: dendropy.Tree) -> None:
    """Checks that every internal node of a tree has two children, naming those that do not."""
    nodes = [node for node in tree.preorder_internal_node_iter() if len(node.child_nodes()) != 2]
    if not nodes:
        return
    problems = []
    for node in nodes:
        children = node.child_nodes()
        problem = f"{len(children)} branch{'es' if len(children) != 1 else ''} to "
        problem += fmt_species_list(first_tip(child) for child in children)
        if node is tree.seed_node:
            problem += " (at the root: if this tree is unrooted it must be rooted first)"
        problems.append(problem)
    msg = (
        f"the tree is not binary; {len(nodes)} node{'s have' if len(nodes) != 1 else ' has'} "
        f"other than two branches below {'them' if len(nodes) != 1 else 'it'}:\n{fmt_bullets(problems)}"
    )
    raise click.BadParameter(msg, ctx, param)


def check_branch_lengths(ctx: click.Context | None, param: click.Parameter | None, tree: dendropy.Tree) -> None:
    """Checks that every branch of a tree below its root has a length."""
    nodes = [node for node in tree.preorder_node_iter() if node is not tree.seed_node and node.edge.length is None]
    if not nodes:
        return
    example = next((node for node in nodes if node.is_leaf()), nodes[0])
    if example.is_leaf():
        where = example.taxon.label
    else:
        where = "the ancestor of " + fmt_species_list(first_tip(child) for child in example.child_nodes())
    verb = "has" if len(nodes) == 1 else "have"
    msg = (
        f"{len(nodes)} of {len(tree.nodes()) - 1} branches {verb} no length, such as the one to {where}; "
        "TACT needs a dated tree"
    )
    raise click.BadParameter(msg, ctx, param)


def check_ultrametric(
    ctx: click.Context | None, param: click.Parameter | None, tree: dendropy.Tree, precision: float
) -> None:
    """Checks that all tips of a tree are the same distance from the root, within a relative `precision`."""
    ultra, ((short_tip, short), (long_tip, long)) = is_ultrametric(tree, precision)
    if ultra:
        return
    # is_ultrametric has set the root distances
    tips = tree.leaf_nodes()
    n_short = sum(not math.isclose(tip.root_distance, long, rel_tol=precision) for tip in tips)
    msg = (
        f"the tree is not ultrametric: tips are between {short} ({short_tip}) and {long} ({long_tip}) "
        f"from the root, and {n_short} of {len(tips)} tips are closer to it than {long_tip}\n"
        "Increase --ultrametricity-precision or use phytools::force.ultrametric in R."
    )
    raise click.BadParameter(msg, ctx, param)


class BackboneCommand(click.Command):
    """Helper class to validate a Click Command that contains a backbone tree.

    At a minimum, the Command must contain a `backbone` parameter, which is validated by `validate_newick`
    and checked to ensure it is a binary tree.

    If the command also contains a `taxonomy` parameter, representing a taxonomic phylogeny,
    every tip of `backbone` must match a species in it, ignoring spaces, underscores and case,
    and the backbone is moved onto the taxonomy's TaxonNamespace. An optional `outgroups`
    parameter may add other taxa not in the `taxonomy`.

    If the command also contains an `ultrametricity_precision` parameter, the
    ultrametricity of the `backbone` is also checked.
    """

    def validate_backbone_variables(self, ctx: click.Context, params: dict[str, Any]) -> dict[str, Any]:
        """Validates variables related to the backbone and taxonomy files."""
        param = next(p for p in self.params if p.name == "backbone")
        backbone = validate_newick(ctx, param, params["backbone"])
        if "taxonomy" in params:
            tn = params["taxonomy"].taxon_namespace
            known = {name_key(taxon.label) for taxon in tn}
            tn.new_taxa(name for name in params.get("outgroups") or [] if name_key(name) not in known)
            tn.is_mutable = False
            missing, respelled = adopt_names(backbone, tn)
            if missing:
                msg = (
                    f"{len(missing)} tip{'s are' if len(missing) != 1 else ' is'} not in the taxonomy:\n"
                    f"{fmt_bullets(missing)}\n"
                    f"Add {'them' if len(missing) != 1 else 'it'} to the taxonomy, "
                    "or list outgroups that should stay out of it with --outgroups."
                )
                raise click.BadParameter(msg, ctx, param)
            if respelled:
                click.echo("Note: these backbone tips are spelled differently in the taxonomy, which TACT uses:")
                click.echo(fmt_bullets([f"{name} => {taxon}" for name, taxon in respelled.items()]))

        check_binary(ctx, param, backbone)
        check_branch_lengths(ctx, param, backbone)
        update_tree_view(backbone)

        if "ultrametricity_precision" in params:
            check_ultrametric(ctx, param, backbone, params["ultrametricity_precision"])

        params["backbone"] = backbone
        return params

    def make_context(
        self, info_name: str | None, args: list[str], parent: click.Context | None = None, **extra: object
    ) -> click.Context:
        """Set up the proper Click context for a command handler."""
        ctx = super().make_context(info_name, args, parent, **extra)
        ctx.params = self.validate_backbone_variables(ctx, ctx.params)
        return ctx
