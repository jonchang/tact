import csv

import dendropy
import click


def wc(filename):
    """
    Get number of lines in a file. This should be relatively cheap,
    and it might put the file in the I/O cache.
    """
    with open(filename, "rU") as rfile:
        for ii, _ in enumerate(rfile):
            pass
        return ii + 1


def build_taxonomic_tree(filename):
    """
    Builds a taxonomic tree given a filename. Last column is assumed to
    be a species name and the file must be sorted. All ranks must nest
    completely within the next highest rank.
    """
    lines = wc(filename)
    with open(filename, "rU") as rfile:
        reader = csv.reader(rfile)

        rank_names = next(reader)
        rank_names.pop()  # assume last column is species name

        tree = dendropy.Tree()
        node = tree.seed_node
        tn = tree.taxon_namespace
        tn.is_mutable = True

        row = next(reader)
        for col in row:
            if col:
                node = node.new_child(taxon=tn.new_taxon(col))

        to_add = list()
        stack = row
        known_nodes = dict()
        with click.progressbar(enumerate(reader), label="Generating taxonomy", length=lines) as rf:
            for idx, row in rf:
                for prev, cur in zip(reversed(stack), reversed(row)):
                    if prev == cur and prev != "" and cur != "":
                        break
                    else:
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
        return tree


@click.command()
@click.argument("taxonomy", type=click.Path(exists=True, dir_okay=False, readable=True))
@click.option("--output", help="name of the output taxonomic tree", required=True, type=click.Path(writable=True))
@click.option("--schema", help="format of the output taxonomic tree", default="newick", type=click.Choice(["newick", "nexus", "nexml"]))
def main(taxonomy, output, schema):
    """Generates a taxonomic tree from TAXONOMY.

    TAXONOMY is formatted as a CSV file where each column is a taxonomic
    rank (from most inclusive to least inclusive, with the last column as the
    species name) and each row is a separate species. The CSV must be sorted.
    Each rank must be named (i.e., there should be no empty cells in the
    spreadsheet) and unique.

    This script makes **many** assumptions about its input for speed. Check
    the example taxonomy in the examples/ folder for guidance.
    """
    taxonomy = build_taxonomic_tree(taxonomy)
    taxonomy.write_to_path(output, schema=schema)
    click.echo("Output written to: %s" % click.format_filename(output))
