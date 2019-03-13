# TACT - Taxonomy addition for complete trees

[![Build Status](https://travis-ci.com/jonchang/tact.svg?token=CAAYReeKsDcnZM7jk2wY&branch=master)](https://travis-ci.com/jonchang/tact)

Adds tips to a backbone phylogeny using taxonomy simulated with birth-death models


# Installation

If you don't use `pipsi`, you're missing out.
Here are [installation instructions](https://github.com/mitsuhiko/pipsi#readme).

Simply run:

    $ pipsi install .


# Example

Files used are in the [examples](https://github.com/jonchang/tact/tree/master/examples) folder.

Build a taxonomic tree using the provided CSV file. Run `tact_build_taxonomic_tree --help` to see the required format for this file.

```console
$ tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
Output written to: Carangaria.taxonomy.tre
```

`Carangaria.taxonomy.tre` now contains a Newick phylogeny with many polytomies and named nodes indicating relevant taxonomic ranks. Now run the TACT stochastic polytomy resolver algorithm in conjunction with the backbone phylogeny `Caragaria.tre`.

```console
$ tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted --verbose --verbose
Rates  [####################################]  226/226
TACT  [####################################]  642/642  Carangaria
```

There will be several files created with the prefix `Carangaria.tacted`. These include `newick.tre` and `nexus.tre` (your primary output in the form of Newick and NEXUS format phylogenies), `rates.csv` (estimated diversification rates on the backbone phylogeny), and `log.txt` (extremely verbose output on what TACT is doing and why).

You should check the TACT results now for any issues:

```console
$ tact_check_results Carangaria.tacted.newick.tre --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre > checkresults.csv
```

Open up `checkresults.csv` in your favorite spreadsheet viewer and check the `warnings` column for any issues.

# Citation

The manuscript for TACT is currently in review.

TACT owes its existence to much foundational work in the area of stochastic polytomy resolution, namely PASTIS and CorSiM.

* Thomas, G. H., Hartmann, K., Jetz, W., Joy, J. B., Mimoto, A., & Mooers, A. O. (2013). PASTIS: an R package to facilitate phylogenetic assembly with soft taxonomic inferences. Methods in Ecology and Evolution, 4(11), 1011–1017. doi:[10.1111/2041-210x.12117](https://doi.org/10.1111/2041-210X.12117)

* Cusimano, N., Stadler, T., & Renner, S. S. (2012). A New Method for Handling Missing Species in Diversification Analysis Applicable to Randomly or Nonrandomly Sampled Phylogenies. Systematic Biology, 61(5), 785–792. doi:[10.1093/sysbio/sys031](https://doi.org/10.1093/sysbio/sys031)
