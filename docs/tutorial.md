# Tutorial: your first TACT analysis

This tutorial walks you through a complete TACT analysis using example data. By the end, you'll know how to use TACT to add missing species to an incomplete phylogeny.

!!! note "What you'll learn"
    - How to prepare your data for TACT
    - How to build a taxonomic tree from a CSV file
    - How to run TACT to add missing species
    - How to check your results for problems

## Getting the example data

The example files are available in the [examples folder](https://github.com/jonchang/tact/tree/HEAD/examples) on GitHub. We'll use data from the Carangaria group of fishes.

Download the example files:

```bash
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.csv
curl -LO https://raw.githubusercontent.com/jonchang/tact/HEAD/examples/Carangaria.tre
```

These files contain:

- `Carangaria.csv`: A **taxonomy** file listing all species and their taxonomic ranks
- `Carangaria.tre`: A **backbone phylogeny** with only some species included

## Building the taxonomic tree

First, we need to convert the CSV taxonomy file into a phylogenetic tree format that TACT can use.

??? tip "How do I format the CSV file?"
    Run `tact_build_taxonomic_tree --help` to see the required format for the CSV file, or read the [command reference](commands.md) for more details. The `Carangaria.csv` shows what the format could look like:
    
    | series     | order         | suborder               | family     | genus   | genus.species       |
    |------------|---------------|------------------------|------------|---------|---------------------|
    | Carangaria | Carangiformes | Carangiformes suborder | Carangidae | Alectis | Alectis alexandrina |
    | Carangaria | Carangiformes | Carangiformes suborder | Carangidae | Alectis | Alectis ciliaris    |
    | Carangaria | Carangiformes | Carangiformes suborder | Carangidae | Alepes  | Alepes apercna      |
    | Carangaria | Carangiformes | Carangiformes suborder | Carangidae | Alepes  | Alepes djedaba      |
    | Carangaria | Carangiformes | Carangiformes suborder | Carangidae | Atropus | Atropus atropos     |


```console
$ tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
Output written to: Carangaria.taxonomy.tre
```

This creates `Carangaria.taxonomy.tre`, a [Newick format](https://en.wikipedia.org/wiki/Newick_format) phylogenetic tree with many polytomies representing taxonomic relationships. Each taxonomic rank (genus, family, etc.) is represented as a named node in the tree.

??? tip "Can I skip this step?"
    Yes. `tact_add_taxa` also accepts the CSV file directly (`--taxonomy Carangaria.csv`) and builds the same tree in memory, so it gives the same results. Building the tree first is still useful if you want to look at it, or to check the results with `tact_check_results` later, which needs the tree.

## Adding unsampled species

Now we'll run TACT to add the missing species from the taxonomy tree onto the backbone phylogeny. This is the main TACT command:

```console
$ tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted
Rates  [############]  227/227
TACT  [############]  642/642  Carangaria
```

You'll see two progress bars:

- **Rates**: TACT is estimating diversification rates for each taxonomic group
- **TACT**: TACT is adding unsampled species to the tree

If the output is not a terminal (for example, in the log of a cluster job), or if you add `--verbose`, TACT prints a line every 10% of the way instead of drawing bars.

TACT places species at random, so each run gives a different tree. To get the same tree again, for example to share an analysis that others can reproduce, add `--seed` with any whole number, such as `--seed 42`, and keep it with the inputs and the TACT version you used.

### Understanding the output files

TACT creates several files with the prefix `Carangaria.tacted` (the argument we passed to `--output`):

- **`Carangaria.tacted.newick.tre`** and **`Carangaria.tacted.nexus.tre`**: Your complete phylogenies in Newick and NEXUS formats. These are your main results—complete, binary trees ready for downstream analysis.
- **`Carangaria.tacted.rates.csv`**: Estimated diversification rates (birth and death rates) for each taxonomic group. This shows what rates TACT used when placing species.
- **`Carangaria.tacted.placements.csv`**: Where each added species was placed. See [below](#species-placed-in-a-higher-taxon).
- **`Carangaria.tacted.log.txt`**: Log of what TACT did and why. Useful for troubleshooting if something goes wrong. Add `-v` to also see these messages as TACT runs, or `-vv` to include debugging messages.

### Species placed in a higher taxon

When TACT finishes, it prints a summary like this:

```
Added 642 species; 241 (38%) were placed in a higher taxon because the lowest taxon with species in the backbone is not monophyletic there.
Most affected:
  Soleidae: 44 species placed in Pleuronectoidei (1 rank up; broken up by Cynoglossus attenuatus, Cynoglossus bilineatus and 17 others)
  Arnoglossus: 30 species placed in Bothidae (1 rank up; broken up by Asterorhombus cocosensis, Bothus lunatus and 18 others)
  ...
Placement of each added species: Carangaria.tacted.placements.csv
```

TACT adds each missing species within the clade of its genus (or, for a genus with no sampled species, the next rank up). Here, the 44 missing species of Soleidae could not be added within Soleidae, because in the backbone Soleidae's clade also contains *Cynoglossus*. They were added anywhere within Pleuronectoidei instead. `Carangaria.tacted.placements.csv` lists every added species, where it was expected and where it went; see [`tact_add_taxa`](commands.md#tact_add_taxa) for the columns. If many species of a taxon you care about were placed higher, consider fixing the taxonomy (for example, by splitting a non-monophyletic genus, or creating a new *sensu lato* clade) before relying on the results.

## Checking your results

You should always check your TACT results for any issues:

```console
$ tact_check_results Carangaria.tacted.newick.tre --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre > checkresults.csv
```

This command validates your results by checking for consistency between:

- The simulated (complete) tree you just created
- The original backbone tree
- The taxonomy tree

### Interpreting the results

Open `checkresults.csv` in a spreadsheet program (Excel, Google Sheets, etc.) and look at the **`warnings`** column. This column will be empty for taxonomic groups that passed all checks. If you see warnings, they might indicate:

- Monophyly issues (groups that should be together but aren't)
- Mismatches between expected and observed species counts
- Other potential problems

For more details on interpreting results and troubleshooting issues, see the [Troubleshooting Guide](troubleshooting.md).

## What's next?

Now that you've completed your first TACT analysis, you can:

- **Use your complete tree**: The `.newick.tre` or `.nexus.tre` files can be used in downstream analyses (like BAMM for detecting diversification rate shifts)
- **Learn more**: Read the [Background](background.md) to understand how TACT works and when to use it
- **Explore commands**: Check the [Command Reference](commands.md) for detailed documentation and advanced options
- **Troubleshoot**: If you encountered issues, see the [Troubleshooting Guide](troubleshooting.md)

!!! warning "Important limitations"
    **Do not use TACT-generated trees for trait evolution analyses.** TACT is designed for diversification rate studies, not trait evolution. See the [Background](background.md#important-limitations) for details.
