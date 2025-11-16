# Introducing TACT

!!! abstract "For more details"
    The best way to understand TACT, why it exists, and how it's used, is in my **[point of view paper in _Systematic Biology_](https://academic.oup.com/sysbio/article/69/3/602/5658637?guestAccessKey=cfe86cce-9d88-4f62-80ce-1f66a8f5ff1d)**. This page provides a brief, accessible overview.

Understanding how species diversify over time—the rates at which new species arise and go extinct—is fundamental to explaining major biodiversity patterns on Earth. These patterns range from the [**latitudinal diversity gradient**](https://en.wikipedia.org/wiki/Latitudinal_diversity_gradient) to temporal patterns of [**mass extinction**](https://en.wikipedia.org/wiki/Extinction_event) and recovery. To study diversification, researchers typically need complete or near-complete [**phylogenies**](https://en.wikipedia.org/wiki/Phylogenetic_tree), evolutionary trees that include all or most species in the group they're studying.

However, in practice, most phylogenies are incomplete, because not all species have been sampled for genetic data. This creates a fundamental problem: how can you accurately estimate [**diversification rates**](https://en.wikipedia.org/wiki/Diversification_rate) when you're missing most of your data?

## The problem: incomplete phylogenies

When you have an incomplete phylogeny, you know that some species are missing, but you don't know where they should be placed on the tree. For example, imagine you're studying a genus that contains 10 species, but your phylogeny only includes 3 of them. Where do the other 7 species belong? What are their relationships to the sampled species?

One approach is to add a **[polytomy](https://en.wikipedia.org/wiki/Polytomy)** at the genus level, indicating that all the unsampled species belong somewhere within that group. However, most analysis tools require **[binary trees](https://en.wikipedia.org/wiki/Phylogenetic_tree#Bifurcating_versus_multifurcating)**, where each internal node splits into exactly two branches. Without additional information, you cannot resolve this polytomy to create a binary tree for further analysis.

This is where stochastic polytomy resolution comes in.

## What is stochastic polytomy resolution?

**Stochastic polytomy resolution** is a method that generates complete, binary phylogenies from incomplete phylogenies, by:

1. Taking the polytomies in your incomplete phylogeny
2. Using information about where unsampled species might belong (typically from taxonomy)
3. Placing these unsampled species onto the tree to resolve the polytomies
4. Generating multiple possible complete trees

Because the true placement of unsampled species is unknown, this last step is crucial. TACT generates a **distribution of possible trees**, with each tree representing one possible way the unsampled species could be placed. You can then account for this uncertainty in your downstream analyses.

This approach is particularly powerful because it allows you to use tools designed for complete phylogenies, even when your original data are incomplete.

## Why use TACT?

When diversification rates vary across a phylogeny, incomplete sampling makes it much harder to detect these shifts in diversification rates. The standard "[sampling fraction](#the-sampling-fraction-method)" method, while unbiased, has **reduced statistical power** when sampling is sparse. This means you might miss important rate shifts that actually exist in your data.

??? warning "Learn more about the statistical power problem"
    When a phylogeny is incomplete and diversification rates vary, your ability to detect rate variation declines dramatically. The sampling fraction method produces "flatter" [likelihood](https://en.wikipedia.org/wiki/Likelihood_function) surfaces with incomplete data, making it difficult to confidently identify where rates have changed. This is a fundamental limitation that affects many common analyses.

TACT addresses this by generating complete phylogenies through stochastic polytomy resolution. By analyzing multiple realizations of where unsampled species might be placed, you can more confidently detect diversification rate variation when it exists.

## When should you use TACT?

The choice of method depends on how complete your phylogeny is. Here are our recommendations:

!!! success "Sampling between 1-50%"
    **This is where TACT is most valuable.** When your phylogeny is moderately incomplete (between 1% and 50% sampling), TACT can significantly improve your ability to detect diversification rate shifts that would be missed using the sampling fraction method alone.

??? tip "Sampling above 50%"
    If your phylogeny includes more than 50% of the species in your clade of interest, the "sampling fraction" method (used in BAMM, BiSSE, RPANDA) is usually sufficient. TACT may not provide significant benefits in this case, and the additional complexity may not be worth it.

??? warning "Sampling below 1%"
    For extremely incomplete phylogenies (less than 1% sampling), the "taxonomic" method (used in MEDUSA, APE) may be more appropriate.

## Alternatives to stochastic polytomy resolution

Two other common methods exist to account for incomplete sampling on phylogenies. Understanding these alternatives helps clarify when TACT is the best choice.

### The sampling fraction method

The "sampling fraction" method accounts for incomplete sampling by including a parameter that represents what fraction of species in a clade are included in the phylogeny. This method is used in tools such as [BAMM](http://bamm-project.org), [BiSSE](https://cran.r-project.org/package=diversitree), and [RPANDA](https://cran.r-project.org/package=RPANDA).

**When it works well:**

- Sampling is relatively complete (above 50%)
- You can accurately estimate sampling fractions for different parts of the tree
- Diversification rates are relatively constant

**Limitations:**

- Loses statistical power when sampling is sparse (below 50%)
- Can be difficult to specify accurate sampling fractions when sampling is taxonomically biased
- Estimating trait-specific diversification rates becomes more challenging with lower sampling ([FitzJohn et al. 2009](https://doi.org/10.1093/sysbio/syp067))
- Specifying clade-specific sampling fractions can lead to improper likelihoods in some cases

### The taxonomic method

The "taxonomic" method uses phylogenies where terminal tips represent exemplar taxa (e.g., a single species representing an entire genus), with species richness assigned to each tip. The stem age of the exemplar is then used to compute diversification rates. This method is used in tools such as [MEDUSA](https://cran.r-project.org/package=geiger) and [APE](https://cran.r-project.org/package=ape).

**When it works well:**

- Extremely incomplete phylogenies (below 1% sampling)
- You have good taxonomic information but poor phylogenetic sampling

**Limitations:**

- Extremely sensitive to exactly which node is selected for estimation ([May and Moore 2016](https://doi.org/10.1093/sysbio/syw026))
- Requires discarding information about branching patterns that occur below the level of exemplar taxa
- Cannot estimate parameters or rate shifts below the level of exemplar representation
- Has been shown to be nonidentifiable in some cases ([Rabosky and Benson 2021](https://doi.org/10.1038/s41467-021-23307-5))

### How TACT compares

TACT and other stochastic polytomy resolution methods address a key weakness of both approaches: **when diversification rates vary and sampling is incomplete, your ability to detect rate variation is greatly reduced**. By generating complete phylogenies through stochastic placement of unsampled taxa, TACT can significantly improve statistical power to detect diversification rate shifts.

The main advantage of TACT over other stochastic polytomy resolvers (like PASTIS or CorSiM) is that TACT estimates diversification rates locally for each taxonomic group, rather than assuming a single rate across the entire tree. This allows for more accurate placement of unsampled taxa when rates vary across the phylogeny.

Ready to try TACT? Check out the **[Tutorial](tutorial.md)** for a step-by-step walkthrough with example data, or see the **[Command Reference](commands.md)** for detailed documentation of all TACT commands.

## Important Limitations

!!! danger "Not for trait evolution"
    **Do not use TACT-generated phylogenies for analyses of trait evolution.** Stochastic polytomy resolution can eliminate true [phylogenetic signal](https://en.wikipedia.org/wiki/Phylogenetic_signal) in traits and inflate estimates of trait evolution rates. This occurs because the polytomy resolver does not consider potential phylogenetic conservatism of trait states, so any given tip will likely be placed in a location that overdisperses trait values. See [Rabosky (2015)](https://doi.org/10.1111/evo.12817) for the full reasoning.

    If you need to study trait evolution with incomplete phylogenies, consider using `locate.yeti` from [Revell et al. (2015)](https://doi.org/10.1111/evo.12628), but this requires having good trait sampling for your clade of interest.

!!! note "What TACT does not account for"

    - **Topological uncertainty** in the backbone phylogeny (the tree structure is assumed to be correct)
    - **Temporal uncertainty** in [node](https://en.wikipedia.org/wiki/Node_(phylogenetics)) ages (node ages are assumed to be known)
    - **Temporal rate variation** (TACT handles rate variation across lineages, but not through time)

    If you suspect your data may fall into one of these scenarios, other methods might be more appropriate.
