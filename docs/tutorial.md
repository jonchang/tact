Files used are in the [examples](https://github.com/jonchang/tact/tree/HEAD/examples) folder.

```console
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.csv
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.tre
```

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
