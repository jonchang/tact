"""Command-line interface module to construct taxonomic phylogenies."""

import csv
import re
from dataclasses import dataclass, field

import click
import dendropy

from .names import clean_name, name_key

# Single root above the highest rank, so that the tree has one root even if the first column has several names
ROOT = "__TAXONOMIC_ROOT__"

# Problems listed in an error message; the rest are counted
MAX_PROBLEMS = 20

# Characters that cannot appear in an unquoted Newick label
_NEEDS_QUOTES = re.compile(r"[()\[\]':;,_\"\t\r\n]")


class TaxonomyError(click.ClickException):
    """A taxonomy CSV that cannot be turned into a tree, with every problem found in it."""

    def __init__(self, problems: list[str]) -> None:
        """Summarize `problems`, a list of messages that start with the line number they refer to."""
        self.problems = problems
        lines = [f"Found {len(problems)} problem{'s' if len(problems) != 1 else ''} in the taxonomy:"]
        lines.extend(f"* {problem}" for problem in problems[:MAX_PROBLEMS])
        if len(problems) > MAX_PROBLEMS:
            lines.append(f"* ... and {len(problems) - MAX_PROBLEMS} more")
        super().__init__("\n".join(lines))


@dataclass
class Clade:
    """A taxon and the taxa directly below it, keyed by name. Species have no `Clade`."""

    children: dict[str, "Clade | None"] = field(default_factory=dict)


@dataclass
class Taxonomy:
    """A taxonomy read from a CSV file.

    Attributes:
        ranks: Names of the rank columns, from most to least inclusive, without the species column.
        root: Taxa of the highest rank.
        renamed: Names used at more than one rank, mapped from `(name, rank)` to the name given at that rank.
        filled: Number of empty cells filled with placeholder taxa.
        respelled: Rank names spelled more than one way, each spelling with the first line it is on, mapped from
            the spelling used.
    """

    ranks: list[str]
    root: Clade
    renamed: dict[tuple[str, str], str]
    filled: int
    respelled: dict[str, list[tuple[str, int]]] = field(default_factory=dict)

    def to_newick(self) -> str:
        """Write the taxonomy as a Newick tree with labelled internal nodes, children in alphabetical order."""

        def walk(name: str, clade: Clade | None) -> str:
            if clade is None:
                return quote(name)
            children = ",".join(walk(child, clade.children[child]) for child in sorted(clade.children))
            return f"({children}){quote(name)}"

        return f"({walk(ROOT, self.root)});"


def quote(label: str) -> str:
    """Format `label` for Newick the way DendroPy reads it back: spaces as underscores, quoted if needed."""
    if _NEEDS_QUOTES.search(label):
        return "'" + label.replace("'", "''") + "'"
    return label.replace(" ", "_")


def _fill_empty(ranks: list[str], species: str, rank_names: list[str]) -> int:
    """Name each empty rank after the nearest named taxon below it and the rank, e.g. `Gobius family`."""
    names = [*ranks, species]
    filled = 0
    for i, name in enumerate(ranks):
        if not name:
            below = next(n for n in names[i + 1 :] if n)
            ranks[i] = clean_name(f"{below} {rank_names[i]}")
            filled += 1
    return filled


def read_taxonomy(path: str, fill_empty: bool = False) -> Taxonomy:
    """Read a taxonomy from a CSV file with one row per species and one column per rank.

    The columns go from the most to the least inclusive rank, and the last column is the species name.
    Rows can be in any order. Every problem in the file is collected before raising.

    * Names that differ only in spaces, underscores or case are the same name, spelled as on its first line
      with single spaces. A species spelled more than one way is an error.
    * A name used at more than one rank (such as a genus and its nominal subgenus) is renamed at the lower
      ranks by appending the rank, as `name rank`. Species names are never renamed.
    * The same name at one rank under two different parents is an error.
    * An empty cell is an error unless `fill_empty` is set, in which case it becomes a placeholder taxon
      named after the nearest named taxon below it.
    * A species listed twice is an error unless both rows are the same.

    Args:
        path: CSV file to read.
        fill_empty: Fill empty rank cells with placeholder taxa instead of reporting them.

    Returns:
        The taxonomy.

    Raises:
        TaxonomyError: If the file has any of the problems above, or rows with the wrong number of cells.
    """
    problems: list[str] = []
    # utf-8-sig drops the byte order mark that spreadsheet programs often add
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    if not rows:
        raise TaxonomyError(["line 1: the file is empty"])

    header = [cell.strip() for cell in rows[0]]
    if len(header) < 2:
        problems.append("line 1: there must be at least one rank column and a species column")
    duplicated = sorted({name for name in header if header.count(name) > 1})
    if duplicated:
        problems.append(f"line 1: duplicated column names: {', '.join(duplicated)}")
    if problems:
        raise TaxonomyError(problems)
    rank_names = header[:-1]
    n_ranks = len(rank_names)

    # Check each row, and note the first line of each spelling of each name
    raw_records: list[tuple[int, list[str], str]] = []
    spellings: dict[str, dict[str, int]] = {}
    species_keys: set[str] = set()
    for lineno, row in enumerate(rows[1:], start=2):
        cells = [cell.strip() for cell in row]
        if not any(cells):
            continue
        if len(cells) != len(header):
            problems.append(f"line {lineno}: expected {len(header)} cells, found {len(cells)}")
            continue
        *ranks, species = cells
        if not species:
            problems.append(f"line {lineno}: no species name")
            continue
        if not all(ranks) and not fill_empty:
            empty = ", ".join(rank_names[i] for i, name in enumerate(ranks) if not name)
            problems.append(
                f"line {lineno}: empty {empty} for {clean_name(species)} (use --fill-empty to fill it in)"
            )
            continue
        for name in cells:
            if name:
                spellings.setdefault(name_key(name), {}).setdefault(name, lineno)
        species_keys.add(name_key(species))
        raw_records.append((lineno, ranks, species))

    # Spell each name as on its first line, and report species spelled more than one way
    respelled: dict[str, list[tuple[str, int]]] = {}
    spelling: dict[str, str] = {}
    for key, found in spellings.items():
        first = next(iter(found))
        spelling[key] = clean_name(first)
        if len(found) > 1 and key not in species_keys:
            respelled[spelling[key]] = list(found.items())
    records: list[tuple[int, list[str], str]] = []
    filled = 0
    for lineno, raw_ranks, raw_species in raw_records:
        found = spellings[name_key(raw_species)]
        first = next(iter(found))
        if raw_species != first:
            if found[raw_species] == lineno:
                problems.append(
                    f"line {lineno}: {raw_species} is spelled {first} on line {found[first]}; "
                    "spell each species the same way on every line"
                )
            continue
        ranks = [spelling[name_key(name)] if name else "" for name in raw_ranks]
        species = spelling[name_key(raw_species)]
        if not all(ranks):
            filled += _fill_empty(ranks, species, rank_names)
        records.append((lineno, ranks, species))

    # Rename names used at more than one rank everywhere except their highest rank, or everywhere if a species
    # has that name, so that which occurrence is renamed does not depend on the order of the rows
    levels: dict[str, set[int]] = {}
    for _, ranks, species in records:
        for i, name in enumerate(ranks):
            levels.setdefault(name, set()).add(i)
        levels.setdefault(species, set()).add(n_ranks)
    renamed: dict[tuple[str, str], str] = {}
    for name, found_at in levels.items():
        if len(found_at) > 1:
            keep = n_ranks if n_ranks in found_at else min(found_at)
            for i in found_at - {keep}:
                renamed[name, rank_names[i]] = clean_name(f"{name} {rank_names[i]}")

    # Insert each row from the highest rank down
    root = Clade()
    # parent and first line of each name at each rank
    parents: list[dict[str, tuple[str, int]]] = [{} for _ in range(n_ranks)]
    reported: set[tuple[int, str, str]] = set()
    species_rows: dict[str, tuple[tuple[str, ...], int]] = {}
    for lineno, ranks, species in records:
        ranks = [renamed.get((name, rank_names[i]), name) for i, name in enumerate(ranks)]
        clade = root
        parent = ROOT
        conflict = False
        for i, name in enumerate(ranks):
            first_parent, first_line = parents[i].setdefault(name, (parent, lineno))
            if first_parent != parent:
                if (i, name, parent) not in reported:
                    reported.add((i, name, parent))
                    problems.append(
                        f"line {lineno}: {rank_names[i]} {name} is in {parent}, "
                        f"but in {first_parent} on line {first_line}"
                    )
                conflict = True
                break
            child = clade.children.setdefault(name, Clade())
            assert child is not None
            clade = child
            parent = name
        if conflict:
            continue
        seen = species_rows.get(species)
        if seen is not None:
            if seen[0] != tuple(ranks):
                problems.append(f"line {lineno}: {species} is also listed on line {seen[1]}, under different ranks")
            continue
        species_rows[species] = (tuple(ranks), lineno)
        clade.children[species] = None

    if problems:
        raise TaxonomyError(problems)
    return Taxonomy(rank_names, root, renamed, filled, respelled)


def echo_changes(taxonomy: Taxonomy) -> None:
    """Tell the user which names were respelled or renamed, and how many empty cells were filled."""
    if taxonomy.respelled:
        click.echo("Note: these names are spelled more than one way, and each is read as the same taxon:")
        for name, found in sorted(taxonomy.respelled.items()):
            spelled = ", ".join(f"{spelling} (line {lineno})" for spelling, lineno in found)
            click.echo(f"{name}: {spelled}")
    if taxonomy.filled:
        click.echo(f"Filled {taxonomy.filled} empty cells with placeholder taxa.")
    if taxonomy.renamed:
        click.echo("Note: several rank names were adjusted to ensure uniqueness. These are:")
        for (name, _), new in sorted(taxonomy.renamed.items()):
            click.echo(f"{name} => {new}")


def build_taxonomic_tree(filename: str, fill_empty: bool = False) -> dendropy.Tree:
    """Build a taxonomic tree from a CSV file. See `read_taxonomy` for the format and checks."""
    newick = read_taxonomy(filename, fill_empty).to_newick()
    return dendropy.Tree.get(data=newick, schema="newick")


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
@click.option(
    "--fill-empty",
    is_flag=True,
    help="fill empty rank cells with a placeholder taxon named after the taxon below it, e.g. Gobius family",
)
def main(taxonomy: str, output: str, schema: str, fill_empty: bool) -> None:
    r"""Generates a taxonomic tree from TAXONOMY.

    TAXONOMY is a comma-separated values file:

    \b
    * Each row represents a single species, in any order.
    * Each column represents a taxonomic rank, from most inclusive to least
      inclusive, with the last column as the species name.
      - e.g. Family,Genus,Species
    * The first row names the ranks; each name must be unique.
    * Every cell must be filled in, unless --fill-empty is given.
    * A name can only appear under one parent at each rank.
      - NO: Gobiidae,Gobius,Gobius niger and Oxudercidae,Gobius,Gobius roulei
    * Underscores are read as spaces, and names that differ only in spaces,
      underscores or case are the same name.

    TACT lists every problem it finds in the file before stopping.
    """
    parsed = read_taxonomy(taxonomy, fill_empty)
    newick = parsed.to_newick()
    if schema == "newick":
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(newick + "\n")
    else:
        dendropy.Tree.get(data=newick, schema="newick").write(path=output, schema=schema)
    echo_changes(parsed)
    click.echo(f"Output written to: {click.format_filename(output)}")
