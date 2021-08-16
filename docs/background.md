# Introducing TACT

The best way to understand TACT, why it exists, and how it's used, is in my [point of view paper in _Systematic Biology_](https://academic.oup.com/sysbio/article/69/3/602/5658637?guestAccessKey=cfe86cce-9d88-4f62-80ce-1f66a8f5ff1d), but I cover these topics briefly on this page.

To analyze diversification over time, using complete phylogenies often leads to the best results. However, it's rare that you have a phylogeny inferred from genetic data for all species in your focal clade of interest.

Given an incomplete phylogeny, you can add polytomies where unsampled taxa might belong (for example, at the crown node of the genus or family). You now have a complete phylogeny with polytomies, and it's likely that whatever downstream analysis tool you're using will be unable to make sense of your phylogeny. Without additional information, you cannot resolve this phylogeny to create a bifurcating tree for further analysis.

Stochastic polytomy resolution can be used to generate completely sampled phylogenies from incomplete phylogenies. It takes the (conceptual) multifurcations in your phylogeny, and resolves them _stochastically_ to generate a pseudoposterior distribution of potential resolutions of a complete phylogeny. This distribution of phylogenies can then be analyzed using other tools that estimate diversification rates.

## Alternatives to stochastic polytomy resolution

Two common methods exist to account for incomplete sampling on phylogenies. The first is what I call the "sampling fraction" method, where skeletal "backbone" trees have unplaced missing taxa assigned to clades on that phylogeny. This is used in tools such as [BAMM](http://bamm-project.org), [BiSSE](https://cran.r-project.org/package=diversitree), and [RPANDA](https://cran.r-project.org/package=RPANDA). The second I call the "taxonomic" method, where terminal exemplar lineages have some amount of species richness assigned to them, and the stem age is used to compute the diversification rate for that lineage. This is used in tools such as [MEDUSA](https://cran.r-project.org/package=geiger) and [APE](https://cran.r-project.org/package=ape).

Both of these approaches suffer from a specific weakness [described in our paper](https://academic.oup.com/sysbio/article/69/3/602/5658637?guestAccessKey=cfe86cce-9d88-4f62-80ce-1f66a8f5ff1d) that TACT and other stochastic polytomy resolution methods can address. When estimating diversification rates on phylogenies, if a phylogeny is incomplete and diversification rates vary, your ability to reject a hypothesis of constant-rate diversification declines precipitously, unless the phylogeny has been augmented with stochastic polytomy resolution.

Other issues with the "sampling fraction" approach have been identified by other practitioners. [FitzJohn et al. (2009)](https://doi.org/10.1093/sysbio/syp067) finds that estimating trait-specific diversification rate becomes more challenging with lower sampling. And an [argument from Beaulieu](https://cran.r-project.org/web/packages/hisse/vignettes/Clade-specific-sampling.pdf) suggests that specifying clade-specific sampling fractions can lead to improper likelihoods.

With respect to the taxonomic approach, the problems are myriad, but they can be extremely sensitive to exactly which node is selected for estimation ([May and Moore 2016](https://doi.org/10.1093/sysbio/syw026)), and may require dropping lots of information about branching patterns that occur below the level of the exemplar taxa chosen. A paper from [Rabosky and Benson (2021)](https://doi.org/10.1038/s41467-021-23307-5) also suggests that this approach is nonidentifiable.

## Problems with stochastic polytomy resolution

Generally speaking, analyses of trait evolution should not use phylogenies generated via stochastic polytomy resolution. See [Rabosky (2015)](https://doi.org/10.1111/evo.12817) for the reasoning. You may be able to get around this using `locate.yeti` from [Revell et al. (2015)](https://doi.org/10.1111/evo.12628), but this requires having decent trait sampling for your clade of interest.
