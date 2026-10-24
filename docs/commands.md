# Command Reference

## `tact_build_taxonomic_tree`

Generates a taxonomic tree from a CSV file. Output contains polytomies representing taxonomic relationships.

**Syntax:**
```bash
tact_build_taxonomic_tree [OPTIONS] TAXONOMY
```

**Required:**

* `TAXONOMY`: Path to CSV file with taxonomic information
* `--output PATH`: Output file name

**Options:**

* `--schema [newick|nexus|nexml]`: Output format (default: `newick`)
* `--fill-empty`: Fill empty rank cells with placeholder taxa instead of stopping (see below)

**CSV format:** The first row names the ranks. Each following row is one species. Columns are taxonomic ranks from most to least inclusive, and the last column is the species name.

**Names:** TACT reads underscores as spaces and treats names case-insensitively, so `Gobius_niger` and `gobius niger` are considered the same species. TACT always writes names with underscores in the files it creates (the trees and the CSV reports), so that they match the tip labels that R and other programs read from the trees.

Some common problems with the taxonomy CSV format:

* **Wrong number of cells:** every row must have as many cells as the header.
* **A name under two parents:** a name can appear at a rank under only one parent. For example, the genus *Gobius* cannot be in both Gobiidae and Oxudercidae. Rename one of them, for example to `Gobius (Oxudercidae)`.
* **Empty cells:** every rank must be filled in. With `--fill-empty`, an empty cell becomes a placeholder taxon named after the nearest named taxon below it and the rank, such as `Gobius family`. Each placeholder holds only that taxon, so unrelated groups are never merged.
* **Conflicting duplicates:** a species listed twice must have the same ranks both times. Identical rows are ignored.
* **A species spelled two ways:** a species must be spelled the same way on every line, for example not as `Gobius niger` on one line and `Gobius_niger` or `gobius niger` on another.

A rank name spelled more than one way, such as `Gobiidae` and `gobiidae`, is read as one taxon spelled as on its first line, and TACT prints a note listing the spellings.

A name used at more than one rank, such as a genus and its nominal subgenus, is renamed at the lower rank by appending the rank name (`Cichla subgenus`), and TACT prints a note listing the renamed taxa. Species names are never renamed.

??? example "Examples of taxonomy CSV rows"

    **Valid examples:**

    * `Cichlidae,Cichla,Cichla temensis`
    * `Actinopterygii,Teleostei,Carangiformes,Carangidae,Alectis,Alectis alexandrina`
    * `Felidae,Panthera,Panthera leo`
    * `Cichlidae,Cichla,Cichla,Cichla temensis`: genus and subgenus share a name, so the subgenus is renamed

    **Invalid examples:**

    * `Cichlidae,,Cichla temensis`: empty rank (genus is missing), unless you use `--fill-empty`
    * `Cichla,Cichlidae,Cichla temensis`: wrong order (genus before family; ranks must go from most to least inclusive)
    * `Cichlidae,Cichla,Cichla temensis,ExtraColumn`: too many columns (species should be last)

**Example**

```bash
tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
```

## `tact_add_taxa`

Main TACT command: adds unsampled taxa to a backbone phylogeny using stochastic polytomy resolution.

**Syntax:**
```bash
tact_add_taxa [OPTIONS]
```

**Required:**

* `--taxonomy FILE`: Taxonomy, either as a CSV file in the [`tact_build_taxonomic_tree`](#tact_build_taxonomic_tree) format or as the Newick tree that command builds from it.
* `--backbone FILENAME`: Backbone phylogeny
* `--output TEXT`: Base name for output files

**Options:**

* `--outgroups TEXT`: Comma-separated list of outgroup taxa to ignore
* `--min-ccp FLOAT` (0-1, default: 0.8): Minimum [crown capture probability](troubleshooting.md#what-is-crown-capture-probability)
* `--yule`: Use Yule (pure-birth) model (extinction = 0)
* `--seed INTEGER`: Seed for the random number generators. Two runs with the same seed, inputs, options and TACT version write identical trees and CSV files. Without it, each run gives a different tree.
* `--ultrametricity-precision FLOAT` (default: 1e-6): [Ultrametricity check precision](troubleshooting.md#tree-is-not-ultrametric)
* `--fill-empty`: For a taxonomy CSV, fill empty rank cells with placeholder taxa, as with [`tact_build_taxonomic_tree --fill-empty`](#tact_build_taxonomic_tree). Not allowed with a taxonomy tree.
* `-v, --verbose`: Show more log messages as TACT runs: `-v` for details of each step, `-vv` for debugging messages

**Output files:**

* `{output}.newick.tre`, `{output}.nexus.tre`: Complete binary ultrametric phylogenies (primary output)
* `{output}.rates.csv`: Estimated diversification rates (columns: `taxon`, `birth`, `death`, `ccp`, `source`)
* `{output}.placements.csv`: Where each added species was placed
* `{output}.log.txt`: Log of what TACT did and why (includes debugging messages with `-vv`)

**Placement report:**

Each added species is expected to go into the lowest taxon containing it that has species in the backbone. This is typically its genus, or for a genus with no sampled species, the next rank up. If that taxon is not monophyletic in the backbone, TACT instead adds the species to the first enclosing taxon that is, anywhere within that taxon's clade. At the end of a run, TACT reports how many species this happened to and the taxa responsible, and writes one row per added species to `{output}.placements.csv`:

* `species`: The added species
* `expected_taxon`: The lowest taxon containing it with species in the backbone
* `placed_taxon`: The taxon whose clade it was added to
* `ranks_up`: How many ranks above `expected_taxon` that is
* `new_clade`: For species in a taxon with no sampled species, the taxon generated as a new clade for them
* `skipped_taxa`: Taxa between `expected_taxon` and `placed_taxon` that were not monophyletic, separated by `;`
* `n_intruders`, `intruders`: Number of backbone tips outside `expected_taxon` that fall inside its clade, and the first five of them

**Example:**
```bash
tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted
```

Or, skipping `tact_build_taxonomic_tree`:
```bash
tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.csv --output Carangaria.tacted
```

## `tact_check_results`

Validates output from `tact_add_taxa` by checking consistency between simulated tree, backbone, and taxonomy.

**Syntax:**
```bash
tact_check_results [OPTIONS] SIMULATED
```

**Required:**

* `SIMULATED`: Path to simulated tree (typically `{output}.newick.tre`)
* `--backbone FILE`: Original backbone phylogeny
* `--taxonomy FILE`: Taxonomic phylogeny

**Options:**

* `--output FILENAME`: Output CSV (default: stdout)
* `--cores INTEGER` (1-8, default: CPU count): Number of parallel cores
* `--chunksize INTEGER`: Nodes per core (auto-calculated if not specified)

**Output CSV columns:** `node`, `taxonomy_tips`, `backbone_tips`, `simulated_tips`, `backbone_monophyletic`, `simulated_monophyletic`, `backbone_birth`, `simulated_birth`, `backbone_death`, `simulated_death`, `warnings`

Check the `warnings` column for issues. Empty = passed all checks. Otherwise, see [Troubleshooting](troubleshooting.md#warnings-from-tact_check_results) for what the warnings mean.

**Example:**
```bash
tact_check_results result.newick.tre --backbone backbone.tre --taxonomy taxonomy.tre --output check.csv
```

!!! note "Next steps"
    For a complete walkthrough with example data, see the [Tutorial](tutorial.md). Or, for more information about TACT's methods and when to use them, see the [Background](background.md).
