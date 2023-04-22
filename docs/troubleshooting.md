# Troubleshooting TACT

## Initial steps

First, ensure that you've [installed the latest version of TACT](install.md). Next, run TACT in extra verbose mode (`-vvv`) and consult the log files. This will often flag any issues around monophyly or other problems that TACT might be struggling with.

## Why is TACT breaking up my monophyletic clades?

There are two things to check here:

1. Is the clade defined in the taxonomy file?
2. Is the clade actually monophyletic?

If the clade isn't defined, TACT won't consider it when restricting placements in monophyletic clades.

If the clade is defined, then check that the clade is actually monophyletic. If you have a clade M comprised of tips A, B, and C, with a most recent common ancestor node X, clade M is monophyletic if and only if the tips descending from node X are exactly A, B, and C.

If you have a single rogue taxon that breaks up an otherwise monophyletic clade, TACT will mark all edges of that clade as being valid placement points for unsampled taxa. This can cause a mostly nice-looking group to turn into a complete mess.

Two solutions are possible here:

1. Drop the rogue taxon and let TACT place it somewhere more appropriate, or
2. Expand the definition of your clade in the taxonomy file to include the rogue taxon.

## Why is TACT so slow?

Using [Carangaria](https://fishtreeoflife.org/taxonomy/series/Carangaria/) (~1000 spp), I measure TACT's speed with [hyperfine](https://github.com/sharkdp/hyperfine):

| Install Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `docker pull jonchang/tact` | 7.880 ± 0.503 | 7.196 | 8.628 | 1.24 ± 0.14 |
| `brew install jonchang/biology/tact` | 6.405 ± 0.335 | 5.904 | 7.133 | 1.01 ± 0.11 |
| `pipx install --python pypy3 tact` | 6.351 ± 0.600 | 5.717 | 7.086 | 1.00 |

Using [Percomorphaceae](https://fishtreeoflife.org/taxonomy/subdivision/Percomorphaceae/) (~17000 spp), I measure TACT's speed with `time(1)`:

* Docker (pypy): ~460s
* pipx (pypy): ~460s
* Homebrew (cpython): ~2500s (5.4x worse)

TACT runs quite quickly for small datasets. However, for large datasets, it can be quite slow, largely due to factors inherent in Python's object model. TACT could be rewritten in a compiled language, or use different data structures, but this would involve a lot of fiddly work that I don't really want to do. Using PyPy is a reasonable compromise for most people's hardware.

On smaller datasets, the overhead of instantiating a Docker container exceeds any speedup to be gained via PyPy. For larger datasets, the speedup to be gained via PyPy vastly outweighs the easier install of the Homebrew method.

On M1 Macbooks, TACT experiences a significant performance penalty (~30-50x) when running under emulation using Docker. Using the `pipx` install method with Python 3.11 is recommended.

## What is crown capture probability?

Crown capture probability is the probability that a sample of `k` taxa from a clade of `n` total taxa includes a root node, under a Yule process. See [Sanderson (1996)](https://doi.org/10.1093/sysbio/45.2.168) for the full manuscript describing this equation.

This probability is used in two key places in TACT as a cutoff, and can be tweaked using the `--min-ccp` parameter to `tact_add_taxa`.

First, TACT will refuse to estimate diversification rates for taxa where the crown capture probability falls below the cutoff. For example, if the known diversity of a genus is 10 species, but the backbone tree only has 3 species, the probability that you've successfully sampled the crown node is approximately 61%. TACT will instead use the next oldest taxonomy node, such as a family or tribe (if defined), to estimate the diversification rate for this group.

Second, TACT will permit stem placement of unsampled taxa if the crown capture probability of the sampled taxon falls below the cutoff. For example, if the known diversity of a genus is 10 species, but the backbone tree only has 3 species, TACT will consider the stem of the group a valid graft point for the remaining 7 species.

The default cutoff of 0.8 was chosen since it balanced between being relatively loose with respect to stem placements, while still giving reasonable results. You can disable this feature altogether by setting it either to 0, meaning the crown is always assumed to be sampled, or to 0.999 (repeating), meaning the crown is always assumed to not be sampled. The latter will tend to produce odd results and is not recommended.

## Taxon has a minimum age constraint A but oldest generated time was B < A

This message appears when TACT is generating an entire monophyletic clade from scratch (e.g. you have a genus that is completely unsampled), but the only valid place to graft this new clade based on the stem age of the clade is incompatible with valid attachment points in the phylogeny.

Suppose you have a family with three genera: _Canis_ (n=5), _Lycalopex_ (n=6), and _Vulpes_ (n=11). If _Canis_ is completely unsampled, and _Lycalopex_ and _Vulpes_ are well-sampled and reciprocally monophyletic, the only valid placement for the stem of _Canis_ would be on the stem of _Lycalopex_, the stem of _Vulpes_, or (potentially) the stem of the entire canine subfamily. All the more recent (tipwards) graft points are part of these monophyletic grops and are therefore not considered when attaching this new clade.

TACT will generate a single new divergence time with a minimum age constraint to replace the original stem age. This could lead to very short internal branches, or unusual placements of new clades.

## Taxon is fully locked, so attaching to stem

This message appears when TACT is generating an entire monophyletic clade from scratch (e.g. you have a genus that is completely unsampled), but the only valid place to graft this new clade is inside another clade.

Suppose you have a family Caninae with two genera, _Canis_ and _Vulpes_, but you've only sampled species from _Canis_. TACT will generate the entire _Vulpes_ clade from scratch, and attempt to attach it somewhere within Caninae, but all the edges subtending the crown node of that group are already part of _Canis_, and TACT will refuse to break the monophyly of that group.

TACT will generate a single new divergence time constrained so that it will attach somewhere on the stem of the sampled clade. This could lead to very short internal branches, or unusual placements of new clades.

## Tree is not ultrametric

TACT tests for the ultrametricity of your phylogeny prior to doing any work on it. If it falls outside a certain precision threshold (by default, 1e-6, or roughly six digits after the decimal point), it'll refuse to run. This is mostly intended as a smoke test to ensure that you haven't accidentally passed a phylogram or extinct-tip phylogeny to TACT.

Otherwise, if the tree is mostly ultrametric, it'll just tweak the ages a bit to make sure everything lines up. This won't really affect your results, and I go into more exhaustive detail in my [blog post on this subject](https://jonathanchang.org/blog/three-ways-to-check-and-fix-ultrametric-phylogenies/).

## Backbone tree is not binary

This sometimes occurs when you have a basal trifurcation in the phylogeny, representing an unrooted tree. You'll need to resolve the root of the phylogeny somehow to root the tree. TACT will generally assume that a basal trifurcation represents a rooted tree, but this doesn't always happen.

