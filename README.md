# TACT - Taxonomy addition for complete trees

[![PyPI](https://img.shields.io/pypi/v/tact.svg)](https://pypi.org/project/tact/)
[![Build status](https://github.com/jonchang/tact/workflows/Python%20package/badge.svg)](https://github.com/jonchang/tact/actions)

Adds tips to a backbone phylogeny using taxonomy simulated with birth-death models

# Installation

TACT requires Python 3. When possible, we recommend using the PyPy 3 implementation as it can significantly speed up TACT analyses, particularly on large datasets. In addition, TACT depends on the click, DendroPy, NumPy, and SciPy packages.

## Homebrew

Using Homebrew is the recommended way to install TACT. [Install Homebrew on macOS](https://brew.sh) or [Install Homebrew on Linux or Windows 10](https://docs.brew.sh/Homebrew-on-Linux). Once Homebrew has been installed, run

    brew install jonchang/biology/tact

## pipx

If you are unable or unwilling to use Homebrew, the next recommended way to install TACT is via `pipx`. [Install `pipx`](https://pipxproject.github.io/pipx/installation/), then run:

    pipx install tact

If you have PyPy3 installed (via e.g. `brew install pypy3` or [some other source](https://www.pypy.org/download.html)) and feel particularly adventurous, you can try to install a faster PyPy version using:

    pipx install --python pypy3 tact

Note that this will take much longer to install, and the installation will almost certainly fail unless you have the proper compilers set up. If it succeeds though, you should see a rather dramatic improvement in TACT's performance.

## Docker

You can also try using the Docker image if you can't get your Python to cooperate. Install [Docker Desktop](https://www.docker.com/products/docker-desktop) and run the following to download the TACT image:

    docker pull jonchang/tact:latest

Then, run TACT from the container image, giving it access to your current working directory:

    mkdir -p examples
    cd examples
    curl -LO https://raw.githubusercontent.com/jonchang/tact/master/examples/Carangaria.csv
    curl -LO https://raw.githubusercontent.com/jonchang/tact/master/examples/Carangaria.tre
    docker run -it -v "$(pwd)":/workdir -w /workdir jonchang/tact tact_build_taxonomic_tree Carangaria.csv --output Carangaria.taxonomy.tre
    docker run -it -v "$(pwd)":/workdir -w /workdir jonchang/tact tact_add_taxa --backbone Carangaria.tre --taxonomy Carangaria.taxonomy.tre --output Carangaria.tacted

One benefit of using the Docker image is that it automatically uses PyPy3, which can speed up your TACT run, but might also be slower on macOS. Here's a screencast of what this looks like:

[![asciicast](https://asciinema.org/a/347571.svg)](https://asciinema.org/a/347571)

## Other

Other ways of installing TACT, including unpacking the tarball somewhere or directly using `pip`, are neither supported nor recommended.

# Example

Files used are in the [examples](https://github.com/jonchang/tact/tree/master/examples) folder.

```console
curl -LO https://raw.githubusercontent.com/jonchang/tact/master/examples/Carangaria.csv
curl -LO https://raw.githubusercontent.com/jonchang/tact/master/examples/Carangaria.tre
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

# Contributing

Development on TACT uses [`poetry`](https://poetry.eustace.io/). Simply clone the repository and install:

```console
$ git clone https://github.com/jonchang/tact.git
$ cd tact
$ poetry install
```

When releasing a new version of tact, run its tests and bump its revision like so:

```console
$ poetry run pytest # optionally with --script-launch-mode=subprocess
$ poetry version patch # or minor, etc.
$ git commit -p
$ git tag VERSION
$ git push
```

A GitHub Actions workflow will build and publish the new version on PyPI.

# Citation

TACT is described more fully in its manuscript. If you use TACT, please cite:

* Chang, J., Rabosky, D. L., & Alfaro, M. E. (2019). Estimating diversification rates on incompletely-sampled phylogenies: theoretical concerns and practical solutions. Systematic Biology. doi:[10.1093/sysbio/syz081](https://doi.org/10.1093/sysbio/syz081)

TACT owes its existence to much foundational work in the area of stochastic polytomy resolution, namely PASTIS and CorSiM.

* Thomas, G. H., Hartmann, K., Jetz, W., Joy, J. B., Mimoto, A., & Mooers, A. O. (2013). PASTIS: an R package to facilitate phylogenetic assembly with soft taxonomic inferences. Methods in Ecology and Evolution, 4(11), 1011–1017. doi:[10.1111/2041-210x.12117](https://doi.org/10.1111/2041-210X.12117)

* Cusimano, N., Stadler, T., & Renner, S. S. (2012). A New Method for Handling Missing Species in Diversification Analysis Applicable to Randomly or Nonrandomly Sampled Phylogenies. Systematic Biology, 61(5), 785–792. doi:[10.1093/sysbio/sys031](https://doi.org/10.1093/sysbio/sys031)

# Sponsorship

Please consider sponsoring the maintenance of TACT via [GitHub Sponsors](https://github.com/sponsors/jonchang).
