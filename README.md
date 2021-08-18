# TACT: Taxonomic Addition for Complete Trees

[![PyPI](https://img.shields.io/pypi/v/tact.svg)](https://pypi.org/project/tact/)
[![Build status](https://github.com/jonchang/tact/workflows/Python%20package/badge.svg)](https://github.com/jonchang/tact/actions)
[![Docker Hub](https://img.shields.io/docker/pulls/jonchang/tact.svg)](https://hub.docker.com/r/jonchang/tact)

TACT is a Python app for stochastic polytomy resolution. It uses birth-death-sampling estimators across an ultrametric phylogeny to generate branching times for unsampled taxa, using taxonomic information to compatibly place new taxa onto a backbone phylogeny.

## Getting started with TACT

* [Installing](https://tact.jonathanchang.org/install/)
* [Introduction and background](https://tact.jonathanchang.org/background/)
* [Example analysis](https://tact.jonathanchang.org/tutorial/)
* [Troubleshooting problems](https://tact.jonathanchang.org/troubleshooting/)

## Citation

TACT is described more fully in its manuscript. If you use TACT, please cite:

* Chang, J., Rabosky, D. L., & Alfaro, M. E. (2019). Estimating diversification rates on incompletely-sampled phylogenies: theoretical concerns and practical solutions. Systematic Biology. doi:[10.1093/sysbio/syz081](https://doi.org/10.1093/sysbio/syz081)

TACT owes its existence to much foundational work in the area of stochastic polytomy resolution, namely PASTIS and CorSiM.

* Thomas, G. H., Hartmann, K., Jetz, W., Joy, J. B., Mimoto, A., & Mooers, A. O. (2013). PASTIS: an R package to facilitate phylogenetic assembly with soft taxonomic inferences. Methods in Ecology and Evolution, 4(11), 1011–1017. doi:[10.1111/2041-210x.12117](https://doi.org/10.1111/2041-210X.12117)

* Cusimano, N., Stadler, T., & Renner, S. S. (2012). A New Method for Handling Missing Species in Diversification Analysis Applicable to Randomly or Nonrandomly Sampled Phylogenies. Systematic Biology, 61(5), 785–792. doi:[10.1093/sysbio/sys031](https://doi.org/10.1093/sysbio/sys031)

## Sponsorship

Please consider sponsoring the ongoing maintenance of TACT via [GitHub Sponsors](https://github.com/sponsors/jonchang).

Initial development was supported by a National Science Foundation Doctoral Dissertation Improvement Grant (DEB-1601830).
