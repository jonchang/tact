# Change History

## tact 0.5.0

* TACT has a new documentation website, available at [tact.jonathanchang.org](https://tact.jonathanchang.org).
* Adds an experimental command, `tact_add_config`. This uses a configuration-based approach to specify nodes of interest where unsampled species will be placed. This feature is currently undocumented and is expected to have many bugs.
* Adds a `--version` option to most commands.
* Uses a new interval bounds checker to ensure that the union of all possible age constraints on a clade is itself an atomic (single) interval, rather than a disjunction of multiple such intervals.
* Checks for a valid taxonomy tree are moved from `tact_build_taxonomic_tree` to `tact_add_taxa`, ensuring that taxonomic trees generated outside of TACT can still be appropriately validated.
* Drops support for Python 3.7.
* Adds support for Python 3.11.
* Updates NumPy to 1.24.
* Updates SciPy to 1.10.
* Updates DendroPy to 4.6.
* Updates the version of PyPy in the Docker image to use Python 3.9.

## tact 0.4.1

* Extreme age ranges when using the Yule or birth-death models should now cause fewer optimization issues (reported by Alexandre Siqueira).

## tact 0.4.0

* Drops support for Python 3.6.
* `tact_add_taxa` gains `--ultrametricity-precision` to control the precision of ultrametricity checks (reported by Miao Sun, [#230](https://github.com/jonchang/tact/issues/230)).

## tact 0.3.4

* Introduces a new dual-optimizer algorithm, which uses simulated annealing to estimate diversification rates when the standard optimizer fails. This should address optimization problems that occur when estimating parameters on particularly species-rich or species-poor groups.
* Rate estimation is now optimized for cherries (by not estimating them at all).
* Improved reporting of which species in the backbone are breaking desired taxonomic monophyly.
* Full support for Python 3.9.

## tact 0.3.3

* TACT now uses DendroPy 4.5.1.

## tact 0.3.2

* Fixes a numerical precision issue in certain phylogenies with zero length branches (reported by Marcio Pie).
* Logs now have a more standardized format.

## tact 0.3.1

* TACT now provides builds via [Docker Hub](https://hub.docker.com/r/jonchang/tact).

## tact 0.3.0

* `tact_build_taxonomic_tree` now sorts its input on the user's behalf (suggested by Marcio Pie).
* `tact_build_taxonomic_tree` automatically generates unique rank names.
* `tact_build_taxonomic_tree` detects and warns on empty input cells.
* `tact_build_taxonomic_tree` checks that the phylogeny it produces is valid.
* `tact_add_taxa` now has fewer annoying warnings.

## tact 0.2.7

* Fixes some DendroPy messages in `tact_add_taxa` that were erroneously passed to the user.

## tact 0.2.6

* Internal automation improvements.

## tact 0.2.5

* Fixes a rare optimization bug when using the Yule model
* Fixes a rare optimization bug when analysing particularly small phylogenies.

## tact 0.2.4

* Internal automation improvements.

## tact 0.2.3

* `tact_add_taxa` now correctly restores terminal settings when quitting (reported by Joseph W. Brown, [#101](https://github.com/jonchang/tact/issues/101)).
* `tact_add_taxa` now assumes in more places that its input trees are rooted.

## tact 0.2.2

* `tact_add_taxa` no longer emits rooting annotations for Newick-format phylogenies (Joseph W. Brown, [#98](https://github.com/jonchang/tact/pull/98)).

* [TACT's manuscript was published!](https://doi.org/10.1093/sysbio/syz081)

## tact 0.2.1

* Updates to TACT's unit tests and dependencies.

## tact 0.2.0

* `tact_add_taxa` gains a `--yule` option for pure-birth rate estimation.
* Fall back to arbitrary-precision math in more circumstances.

## tact 0.1.4

* Update NumPy dependency to 1.17.

## tact 0.1.3

* Migrate to Poetry build system
* Remove poor-performing parallel MRCA rate calculation algorithm

## tact 0.1.2

* Adds example taxonomy and backbone files to the distribution.
* This is the version that was reviewed for [Chang et al. (2019)](https://doi.org/10.1093/sysbio/syz081).

## tact 0.1.1

* Fixes a phylogeny generation bug in `tact_build_taxonomic_tree`.

## tact 0.1.0

* Initial release.
* This is the version of TACT used in [Rabosky et al. (2018)](https://doi.org/10.1038/s41586-018-0273-1).
