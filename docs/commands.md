# Command reference

## `tact_build_taxonomic_tree`

Generates a taxonomic tree from TAXONOMY.

TAXONOMY is a comma-separated values file with several requirements:

* Each row represents a single species.
* Each column represents a taxonomic rank.
* The columns must be arranged from most inclusive to least inclusive, with the last column as the species name.
    * e.g. Family,Genus,Species
* Each rank must be named
* Each rank must be unique
    * OK: Cichlidae,Cichla,Cichla temensis
    * NO: Cichlidae,Cichlidae,Cichla temensis
    * NO: Cichlidae,,Cichla temensis

See the [taxonomy CSV file in the examples/ folder](https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.csv) for a working example.

## `tact_add_taxa`

Adds unsampled taxa from a complete taxonomic phylogeny (generated from `tact_build_taxonomic_tree`) onto a backbone phylogeny.

## `tact_check_results`

Checks the result of the output from `tact_add_taxa`.
