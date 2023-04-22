import collections
import csv

import click
import dendropy

from .validation import validate_tree_node_depths


def fix_file(filename):
    """
    Slurps a file, and does various checks and fixes:
    * Sorts the file
    * Ensures column names are unique
    """
    with open(filename, "r", encoding="utf-8") as rfile:
        lines = rfile.readlines()

    heads = collections.defaultdict(int)
    for item in lines[0].split(","):
        heads[item] += 1
    bad_heads = []
    for key, value in heads.items():
        if value > 1:
            bad_heads.append(key)

    if len(bad_heads) > 0:
        dupe_strs = "\n*  ".join(bad_heads)
        raise click.UsageError(f"CSV headers must have unique names. Duplicated column names:\n*  {dupe_strs}")

    return [lines[0], *sorted(lines[1:])]


def ensure(st, ctx=""):
    "Ensures that a cell is not empty."
    if len(st) == 0:
        if len(ctx) > 0:
            text = f" Offending line:\n{','.join(ctx)}"
        raise click.UsageError("All cells in the CSV must be nonempty." + text)


def mangle_rank(row, names):
    seen = set()
    new = []
    for idx, item in enumerate(row):
        if item in seen:
            if idx >= len(names):
                dupe_str = ",".join(row)
                raise click.UsageError(f"Invalid name (possibly duplicated): {names}\nSeen at: {dupe_str}")
            item = item + "__" + names[idx] + "__"
        seen.add(item)
        new.append(item)
    return new


def build_taxonomic_tree(filename):
    """
    Builds a taxonomic tree given a filename. Last column is assumed to
    be a species name. All ranks must nest completely within the next
    highest rank.
    """
    lines = fix_file(filename)
    reader = csv.reader(lines)

    rank_names = next(reader)
    rank_names.pop()  # assume last column is species name
    rank_names = ["__ROOT__", *rank_names]

    tree = dendropy.Tree()
    node = tree.seed_node
    tn = tree.taxon_namespace
    tn.is_mutable = True

    # Mangle and uniquify the first row
    mangled_ranks = set()
    row = ["__TAXONOMIC_ROOT__", *next(reader)]
    mangled_row = mangle_rank(row, rank_names)
    for orig, new in zip(row, mangled_row):
        if orig != new:
            mangled_ranks.add((orig, new))
    row = mangled_row

    # Initialize the tree structure
    for col in row:
        ensure(col, ctx=row)
        node = node.new_child(taxon=tn.new_taxon(col))

    to_add = []
    stack = row
    known_nodes = {}
    with click.progressbar(enumerate(reader), width=12, label="Generating taxonomy", length=len(lines)) as rf:
        for _, row in rf:
            # Uniquify row names
            row = ["__TAXONOMIC_ROOT__", *row]
            mangled_row = mangle_rank(row, rank_names)
            for orig, new in zip(row, mangled_row):
                if orig != new:
                    mangled_ranks.add((orig, new))
            row = mangled_row

            prev = None
            for prev, cur in zip(reversed(stack), reversed(row)):
                ensure(cur, ctx=row)
                if prev == cur:
                    break
                to_add.append(cur)
            if prev in known_nodes:
                node = known_nodes[prev]
            else:
                node = tree.find_node_with_taxon_label(prev)
                known_nodes[prev] = node
            while to_add:
                nt = to_add.pop()
                if nt:
                    node = node.new_child(taxon=tn.new_taxon(nt))
            stack = row
    tn.is_mutable = False
    if len(mangled_ranks) > 0:
        click.echo("Note: several rank names were adjusted to ensure uniqueness. These are:")
        for orig, new in mangled_ranks:
            click.echo(f"{orig} => {new}")
    return tree


@click.command()
@click.version_option(package_name="tact")
@click.argument("taxonomy", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("--output", help="name of the output taxonomic tree", required=True, type=click.Path(writable=True))
@click.option(
    "--schema",
    help="format of the output taxonomic tree",
    default="newick",
    type=click.Choice(["newick", "nexus", "nexml"]),
)
def main(taxonomy, output, schema):
    """Generates a taxonomic tree from TAXONOMY.

    TAXONOMY is a comma-separated values file with several requirements.

    \b
    * Each row represents a single species.
    * Each column represents a taxonomic rank.
    * The columns must be arranged from most inclusive to least inclusive,
      with the last column as the species name.
      - e.g. Family,Genus,Species
    * Each rank must be named
    * Each rank must be unique
      - OK: Cichlidae,Cichla,Cichla temensis
      - NO: Cichlidae,Cichlidae,Cichla temensis
      - NO: Cichlidae,,Cichla temensis

    This script makes **many** assumptions about its input for speed. Check
    the example taxonomy in the examples/ folder for guidance.
    """
    taxonomy = build_taxonomic_tree(taxonomy)
    msg = validate_tree_node_depths(None, None, taxonomy)
    if msg:
        click.echo(msg)
    taxonomy.write_to_path(output, schema=schema)
    click.echo(f"Output written to: {click.format_filename(output)}")
