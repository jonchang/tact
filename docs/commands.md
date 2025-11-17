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

**CSV format:** Each row = one species; columns = taxonomic ranks (most to least inclusive, last column = species name). Each rank must be named and unique within a row.

??? example "Examples of taxonomy CSV files"

    **Valid examples:**
    
    * `Cichlidae,Cichla,Cichla temensis`
    * `Actinopterygii,Teleostei,Carangiformes,Carangidae,Alectis,Alectis alexandrina`
    * `Felidae,Panthera,Panthera leo`
    
    **Invalid examples:**
    
    * `Cichlidae,Cichlidae,Cichla temensis` - duplicate rank (Cichlidae appears twice)
    * `Cichlidae,,Cichla temensis` - empty rank (genus is missing)
    * `Cichla,Cichlidae,Cichla temensis` - wrong order (genus before family - ranks must be from most inclusive to least inclusive)
    * `Cichlidae,Cichla,Cichla temensis,ExtraColumn` - too many columns (species should be last)

**Example**

```bash
tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
```

!!! warning "Input format"
    This script makes many assumptions about input for speed. See the [example file](https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.csv) for guidance.

## `tact_add_taxa`

Main TACT command: adds unsampled taxa to a backbone phylogeny using stochastic polytomy resolution.

**Syntax:**
```bash
tact_add_taxa [OPTIONS]
```

**Required:**

* `--taxonomy FILENAME`: Taxonomy tree (from `tact_build_taxonomic_tree`)
* `--backbone FILENAME`: Backbone phylogeny
* `--output TEXT`: Base name for output files

**Options:**

* `--outgroups TEXT`: Comma-separated list of outgroup taxa to ignore
* `--min-ccp FLOAT` (0-1, default: 0.8): Minimum [crown capture probability](troubleshooting.md#what-is-crown-capture-probability)
* `--yule`: Use Yule (pure-birth) model (extinction = 0)
* `--ultrametricity-precision FLOAT` (default: 1e-6): [Ultrametricity check precision](troubleshooting.md#tree-is-not-ultrametric)
* `-v, --verbose`: Verbose output (repeatable: `-v`, `-vv`, `-vvv`)

**Output files:**

* `{output}.newick.tre`, `{output}.nexus.tre`: Complete binary ultrametric phylogenies (primary output)
* `{output}.rates.csv`: Estimated diversification rates (columns: `taxon`, `birth`, `death`, `ccp`, `source`)
* `{output}.log.txt`: Detailed log (use `-vv` or `-vvv` for more detail)

**Example:**
```bash
tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted -vv
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

Check the `warnings` column for issues. Empty = passed all checks.

**Example:**
```bash
tact_check_results result.newick.tre --backbone backbone.tre --taxonomy taxonomy.tre --output check.csv
```

!!! note "Next steps"
    For a complete walkthrough with example data, see the [Tutorial](tutorial.md). Or, for more information about TACT's methods and when to use them, see the [Background](background.md).
