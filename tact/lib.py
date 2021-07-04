# -*- coding: utf-8 -*-
from __future__ import division

import collections
import random
import sys
from decimal import Decimal as D
from math import exp
from math import log
from math import isclose

import dendropy
import numpy as np
from scipy.optimize import minimize, dual_annealing

# Raise on overflow
np.seterr(all="raise")


def get_bd(r, a):
    """Converts turnover and relative extinction to birth and death rates."""
    return -r / (a - 1), -a * r / (a - 1)


def get_ra(b, d):
    """Converts birth and death to turnover and relative extinction rates."""
    return (b - d, d / b)


def wrapped_lik_constant(x, sampling, ages):
    """Wrapper for birth-death likelihood to make optimizing more convenient."""
    return lik_constant(get_bd(*x), sampling, ages)


def wrapped_lik_constant_yule(x, sampling, ages):
    """Wrapper for Yule likelihood to make optimizing more convenient."""
    return lik_constant(get_bd(x[0], 0), sampling, ages)


def two_step_optim(func, x0, bounds, args):
    """Tries to optimize function using the fast L-BFGS-B method, and if that fails, use simulated annealing."""
    result = minimize(func, x0=x0, bounds=bounds, args=args, method="L-BFGS-B")
    if result["success"]:
        return result["x"].tolist()

    result = dual_annealing(func, x0=x0, bounds=bounds, args=args)
    if result["success"]:
        return result["x"].tolist()

    raise Exception(f"Optimization failed: {result['message']} (code {result['status']})")


def optim_bd(ages, sampling, min_bound=1e-9):
    """Optimizes birth death using Scipy"""
    if max(ages) < 0.000001:
        init_r = 1e-3
    else:
        # Magallon-Sanderson crown estimator
        init_r = (log((len(ages) + 1) / sampling) - log(2)) / max(ages)
        init_r = max(1e-3, init_r)
    bounds = ((min_bound, 100), (0, 1 - min_bound))
    result = two_step_optim(wrapped_lik_constant, x0=(init_r, min_bound), bounds=bounds, args=(sampling, ages))
    return get_bd(*result)


def optim_yule(ages, sampling, min_bound=1e-9):
    """Optimizes a Yule model using Scipy"""
    if max(ages) < 0.000001:
        init_r = 1e-3
    else:
        # Magallon-Sanderson crown estimator
        init_r = (log((len(ages) + 1) / sampling) - log(2)) / max(ages)
        init_r = max(1e-3, init_r)
    bounds = ((min_bound, 100), (0, 1 - min_bound))
    result = two_step_optim(wrapped_lik_constant_yule, x0=(init_r, 0.0), bounds=bounds, args=(sampling, ages))
    return get_bd(*result)


def get_lik(vec, rho, x):
    l = vec[0]
    m = vec[1]
    root = 1
    lik1 = (root + 1) * np.log(p1(x[0], l, m, rho))
    lik2 = np.sum(np.log(l * p1(x[1:], l, m, rho)))
    lik3 = -(root + 1) * np.log(1 - p0(x[0], l, m, rho))
    return lik1 + lik2 + lik3


def p0_exact(t, l, m, rho):
    t = D(t)
    l = D(l)
    m = D(m)
    rho = D(rho)
    return D(1) - rho * (l - m) / (rho * l + (l * (D(1) - rho) - m) * (-(l - m) * t).exp())


def p0(t, l, m, rho):
    try:
        return 1 - rho * (l - m) / (rho * l + (l * (1 - rho) - m) * exp(-(l - m) * t))
    except FloatingPointError:
        return float(p0_exact(t, l, m, rho))


def p1_exact(t, l, m, rho):
    """Exact version of p1 using Decimal math."""
    t = D(t)
    l = D(l)
    m = D(m)
    rho = D(rho)
    num = rho * (l - m) ** D(2) * (-(l - m) * t).exp()
    denom = (rho * l + (l * (1 - rho) - m) * (-(l - m) * t).exp()) ** D(2)
    return num / denom


def p1_orig(t, l, m, rho):
    """Original version of p1, here for testing and comparison purposes."""
    try:
        num = rho * (l - m) ** 2 * np.exp(-(l - m) * t)
        denom = (rho * l + (l * (1 - rho) - m) * np.exp(-(l - m) * t)) ** 2
        return num / denom
    except (OverflowError, FloatingPointError):
        return float(p1_exact(t, l, m, rho))


def p1(t, l, m, rho):
    """
    Optimized version of p1 using common subexpression elimination and strength reduction from
    exponentiation to multiplication.
    """
    try:
        ert = np.exp(-(l - m) * t, dtype=np.float64)
        num = rho * (l - m) ** 2 * ert
        denom = (rho * l + (l * (1 - rho) - m) * ert) ** 2
        return num / denom
    except (OverflowError, FloatingPointError):
        return float(p1_exact(t, l, m, rho))


def intp1_exact(t, l, m):
    """Exact version of intp1 using Decimal math."""
    l = D(l)
    m = D(m)
    t = D(t)
    num = D(1) - (-(l - m) * t).exp()
    denom = l - m * (-(l - m) * t).exp()
    return num / denom


def intp1(t, l, m):
    try:
        return (1 - exp(-(l - m) * t)) / (l - m * exp(-(l - m) * t))
    except OverflowError:
        return float(intp1_exact(t, l, m))


def lik_constant(vec, rho, t, root=1, survival=1, p1=p1):
    """
    Calculates the likelihood of a constant-rate birth-death process, conditioned
    on the waiting times of a phylogenetic tree and degree of incomplete sampling.

    Based off of the R function TreePar::LikConstant written by Tanja Stadler.

    T. Stadler. On incomplete sampling under birth-death models and connections
    to the sampling-based coalescent. Jour. Theo. Biol. 261: 58-66, 2009.

    Positional arguments:
    vec -- a two element list of birth and death
    rho -- sampling fraction
    t -- vector of waiting times

    Keyword arguments:
    root -- include the root or not? (default: 1)
    survival -- assume survival of the process (default: 1)

    Returns a likelihood. Or FLOAT_MAX.
    """
    try:
        l = vec[0]
        m = vec[1]
        t.sort(reverse=True)
        lik = (root + 1) * log(p1(t[0], l, m, rho))
        for tt in t[1:]:
            lik += log(l) + log(p1(tt, l, m, rho))
        if survival == 1:
            lik -= (root + 1) * log(1 - p0(t[0], l, m, rho))
        return -lik
    except ValueError:
        return sys.float_info.max


def crown_capture_probability(n, k):
    """
    Calculate the probability that a sample of `k` taxa from a clade
    of `n` total taxa includes a root node, under a Yule process.

    This equation is taken from:

    Sanderson, M. J. 1996. How many taxa must be sampled to identify
    the root node of a large clade? Systematic Biology 45:168-173
    """
    if n < k:
        raise Exception(f"n must be greater than or equal to k (n={n}, k={k})")
    if n == 1 and k == 1:
        return 0  # not technically correct but it works for our purposes
    return 1 - 2 * (n - k) / ((n - 1) * (k + 1))


def get_monophyletic_node(tree, species):
    """Returns the node or None that is the MRCA of the `species` in `tree`."""
    mrca = tree.mrca(taxon_labels=species)
    if not mrca:
        return None
    if mrca and species.issuperset(get_tip_labels(mrca)):
        return mrca


def get_birth_death_rates(node, sampfrac, yule=False, include_root=False):
    """
    Estimates the birth and death rates for the subtree descending from
    `node` with sampling fraction `sampfrac`. Optionally restrict to a
    Yule pure-birth model.
    """
    if yule:
        return optim_yule(get_ages(node, include_root), sampfrac)
    else:
        return optim_bd(get_ages(node, include_root), sampfrac)


def get_ages(node, include_root=False):
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    if include_root:
        ages += [node.age]
    return ages


def get_tip_labels(tree_or_node):
    try:
        return set([x.taxon.label for x in tree_or_node.leaf_node_iter()])
    except AttributeError:
        return set([x.taxon.label for x in tree_or_node.leaf_iter()])


def edge_iter(node, filter_fn=None):
    """
    Iterates over the child edge of `node` and all its descendants.
    Can optionally be filtered by `filter_fn`.
    """
    stack = list(node.child_edge_iter())
    while stack:
        edge = stack.pop()
        if filter_fn is None or filter_fn(edge):
            yield edge
        stack.extend(edge.head_node.child_edge_iter())


def get_tree(path, namespace=None):
    """
    Gets a DendroPy tree from a path and precalculate its node ages and bipartition bitmask.
    """
    tree = dendropy.Tree.get_from_path(
        path, schema="newick", taxon_namespace=namespace, rooting="default-rooted"
    )
    update_tree_view(tree)
    return tree


def update_tree_view(tree):
    """
    Mutates a DendroPy tree object with updated node ages and bipartition bitmask. We also
    correct for minor ultrametricity errors.

    Returns a list of tip labels.
    """
    tree.calc_node_ages(is_force_max_age=True)
    tree.update_bipartitions()
    return get_tip_labels(tree)


def is_binary(node):
    """Is the subtree under `node` a fully bifurcating tree?"""
    for x in node.preorder_internal_node_iter():
        if len(x.child_nodes()) != 2:
            return False
    return True


def is_ultrametric(tree, tolerance=1e-6):
    """Is the `tree` ultrametric, within a specified `tolerance`?

    Uses the relative difference between minimum and maximum root-to-tip distances.
    """
    tree.calc_node_root_distances()
    lengths = {}
    for leaf in tree.leaf_node_iter():
        lengths[leaf.taxon.label] = leaf.root_distance
    t_min = min(lengths.items(), key=lambda x: x[1])
    t_max = max(lengths.items(), key=lambda x: x[1])
    return (isclose(t_min[1], t_max[1], rel_tol=tolerance), (t_min, t_max))


def get_short_branches(node):
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge


def compute_node_depths(tree):
    res = dict()
    for leaf in tree.leaf_node_iter():
        cnt = 0
        for anc in leaf.ancestor_iter():
            if anc.label:
                cnt += 1
        res[leaf.taxon.label] = cnt
    return res


def ensure_tree_node_depths(tree):
    node_depths = compute_node_depths(tree)
    stats = collections.defaultdict(int)
    for v in node_depths.values():
        stats[v] += 1
    msg = ""
    if len(stats) > 1:
        msg += "The tips of your taxonomy tree do not have equal numbers of ranked clades in their ancestor chain:\n"
        for k in sorted(stats.keys()):
            msg += f"* {stats[k]} tips have {k} ranked ancestors\n"
    return msg


# TODO: This could probably be optimized
def get_new_times(ages, birth, death, missing, told=None, tyoung=None):
    """
    Simulates new speciation events in an incomplete phylogeny assuming a
    constnat-rate birth-death process.

    Adapted from the R function TreeSim::corsim written by Tanja Stadler.

    N. Cusimano, T. Stadler, S. Renner. A new method for handling missing
    species in diversification analysis applicable to randomly or
    non-randomly sampled phylogenies. Syst. Biol., 61(5): 785-792, 2012.

    Positional arguments:
    ages -- vector of waiting times
    birth -- birth rate
    death -- death rate
    missing -- number of missing taxa to simulate

    Keyword arguments:
    told -- maximum simulated age (default: `max(ages)`)
    tyoung -- minimum simulated age bound (default: `0`)

    Returns a vector of simulated waiting times.
    """
    if told is None:
        told = max(ages)
    if len(ages) > 0:
        if max(ages) > told and abs(max(ages) - told) > sys.float_info.epsilon:
            raise Exception("Zero or negative branch lengths detected in backbone phylogeny")
    if tyoung is None:
        tyoung = 0

    ages.sort(reverse=True)
    times = [x for x in ages if x <= told and x >= tyoung]
    times = [told] + times + [tyoung]
    ranks = range(0, len(times))
    only_new = list()
    while missing > 0:
        if len(ranks) > 2:
            distrranks = list()
            for i in range(1, len(ranks)):
                temp = ranks[i] * (
                    intp1(times[i - 1], birth, death) - intp1(times[i], birth, death)
                )
                distrranks.append(temp)
            try:
                dsum = sum(distrranks)
                distrranks = [x / dsum for x in distrranks]
                for i in range(1, len(distrranks)):
                    distrranks[i] = distrranks[i] + distrranks[i - 1]
                r = random.uniform(0, 1)
                addrank = min([idx for idx, x in enumerate(distrranks) if x > r])
            except ZeroDivisionError:
                addrank = 0
            except ValueError:
                addrank = 0
        else:
            addrank = 0
        r = random.uniform(0, 1)
        const = intp1(times[addrank], birth, death) - intp1(times[addrank + 1], birth, death)
        try:
            temp = intp1(times[addrank + 1], birth, death) / const
        except ZeroDivisionError:
            temp = 0.0
        xnew = (
            1
            / (death - birth)
            * log((1 - (r + temp) * const * birth) / (1 - (r + temp) * const * death))
        )
        only_new.append(xnew)
        missing -= 1
    only_new.sort(reverse=True)
    return only_new
