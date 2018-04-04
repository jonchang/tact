# -*- coding: utf-8 -*-

from __future__ import division

import csv
import sys
import random
from math import log, exp, ceil
from decimal import Decimal as D
import itertools
import subprocess

import dendropy
import numpy as np

# Initialize grid for birthdeath grid search
births = np.linspace(sys.float_info.epsilon, 5, num=100)
deaths = np.linspace(0, 5, num=100)
params = [(x, y) for (x, y) in itertools.product(births, deaths) if x > y]

def get_bd(r, a):
    """Converts turnover and relative extinction to birth and death rates."""
    r = float(r)
    b = r / (1 + float(a))
    d = b - r
    return b, d

def optim_bd_r(ages, sampling):
    """Optimizes birth death using TreePar and R"""
    script = """cat(optim(c({birth}, {death}), function(v, ...) TreePar::LikConstant(v[1], v[2], ...), x = c({ages}), sampling = {sampling}, lower=c(.Machine$double.xmin,0), method = "L-BFGS-B")$par, " dum")"""
    fmt = script.format(birth=1, death=0.02, ages=",".join([str(x) for x in ages]), sampling = sampling)
    output = subprocess.check_output(["Rscript", "--vanilla", "--default-packages=base,stats", "-e", fmt], stderr=subprocess.STDOUT)
    b, d, _ = output.split(None, 2)
    return float(b), float(d)

def optim_bd_grid(ages, sampling):
    """Optimizes birth death using a grid search"""
    res = [lik_constant(x, sampling, ages) for x in params]
    #import pdb; pdb.set_trace()
    return params[np.argmin(res)]

def update_multiplier_freq(q, d=1.1):
    u = np.random.uniform(0, 1, 2)
    l = 2*log(d)
    m = np.exp(l*(u-.5))
    new_q = q * m
    return new_q

def optim_bd(ages, sampling):
    """Optimizes birth death using a cheap MCMC-like algorithm"""
    new_vec = [0, 0]
    ages = np.sort(np.array(ages))[::-1]
    mm = max(ages)
    if mm < 0.0000001:
        vec = [1, 0.02]
    else:
        vec = [len(ages) / mm, 0.02]
    likelihood = get_lik(vec, sampling, ages)
    for x in xrange(100):
        new_vec = update_multiplier_freq(vec)
        new_likelihood = get_lik(new_vec, sampling, ages)
        if (new_likelihood - likelihood) * 100 >= log(np.random.random()):
            likelihood = new_likelihood
            vec = new_vec
    #import pdb; pdb.set_trace()
    return vec

def optim_bd_scipy(ages, sampling):
    """Optimizes birth death using Scipy"""
    from scipy.optimize import minimize
    return minimize(lambda x: lik_constant(x, sampling, ages), [1., 0.02], bounds=((sys.float_info.epsilon, None), (0, None)), method="TNC")["x"].tolist()


def get_lik(vec, rho, x):
    l = vec[0]
    m = vec[1]
    root=1
    lik1= (root + 1) * np.log(p1(x[0], l, m, rho))
    lik2= np.sum(np.log(l * p1(x[1:], l, m, rho)))
    lik3= - (root + 1) * np.log(1 - p0(x[0], l, m, rho))
    return lik1+lik2+lik3

def lik_constant(vec, rho, t, root=1, survival=1):
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
        for i in range(1, len(t)):
            lik += log(l) + log(p1(t[i], l, m, rho))
        if survival == 1:
            lik -= (root + 1) * log(1 - p0(t[0], l, m, rho))
        return -lik
    except ValueError:
        return sys.float_info.max

def p0_exact(t, l, m, rho):
    t = D(t)
    l = D(l)
    m = D(m)
    rho = D(rho)
    return D(1) - rho * (l - m) / (rho * l + (l * (D(1) - rho) - m) * (-(l - m) * t).exp())

def p0(t, l, m, rho):
    try:
        return 1 - rho * (l - m) / (rho * l + (l * (1 - rho) - m) * exp(-(l - m) * t))
    except OverflowError:
        return float(p0_exact(t, l, m, rho))

def p1_exact(t, l, m, rho):
    """Exact version of p1 using Decimal math."""
    t = D(t)
    l = D(l)
    m = D(m)
    rho = D(rho)
    return rho*(l-m)**D(2) * (-(l-m)*t).exp()/(rho*l+(l*(1-rho)-m)*(-(l-m)*t).exp())**D(2)

def p1(t, l, m, rho):
    try:
        return rho*(l-m)**2 * np.exp(-(l-m)*t)/(rho*l+(l*(1-rho)-m)*np.exp(-(l-m)*t))**2
    except OverflowError:
        return float(p1_exact(t, l, m, rho))

def intp1_exact(t, l, m):
    """Exact version of intp1 using Decimal math."""
    l = D(l)
    m = D(m)
    t = D(t)
    return (D(1) - (-(l - m) * t).exp())/(l - m * (-(l - m) * t).exp())

def intp1(t, l, m):
    try:
        return (1 - exp(-(l - m) * t))/(l - m * exp(-(l - m) * t))
    except OverflowError:
        return float(intp1_exact(t, l, m))


def crown_capture_probability(n, k):
    """
    Calculate the probability that a sample of `k` taxa from a clade
    of `n` total taxa includes a root node, under a Yule process.

    This equation is taken from:

    Sanderson, M. J. 1996. How many taxa must be sampled to identify
    the root node of a large clade? Systematic Biology 45:168-173
    """
    if n < k:
        raise Exception("n must be greater than or equal to k (n={}, k={})".format(n, k))
    return 1 - 2 * (n - k) / ((n - 1) * (k + 1))

def get_monophyletic_node(tree, species):
    """Returns the node or None that is the MRCA of the `species` in `tree`."""
    mrca = tree.mrca(taxon_labels=species)
    if not mrca:
        return None
    if mrca and species.issuperset(get_tip_labels(mrca)):
        return mrca

def get_birth_death_for_node(node, sampfrac):
    """
    Estimates the birth and death rates for the subtree descending from
    `node` with sampling fraction `sampfrac`.
    """
    ages = [x.age for x in node.ageorder_iter(include_leaves=False, descending=True)]
    ages += [node.age]
    return optim_bd(ages, sampfrac)

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
    tree = dendropy.Tree.get_from_path(path, schema="newick", taxon_namespace=namespace)
    tree.calc_node_ages()
    tree.encode_bipartitions()
    return tree

def is_binary(node):
    """Is the subtree under `node` a fully bifurcating tree?"""
    for x in node.preorder_internal_node_iter():
        if len(x.child_nodes()) != 2:
            return False
    return True

def get_short_branches(node):
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge

# From TreeSim::corsim
# N. Cusimano, T. Stadler, S. Renner. A new method for handling missing species in diversification analysis applicable to randomly or non-randomly sampled phylogenies. Syst. Biol., 61(5): 785-792, 2012.
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
    assert max(ages) <= told
    if tyoung is None:
        tyoung = 0

    ages.sort(reverse=True)
    n = len(ages) + 1
    times = [x for x in ages if x <= told and x >= tyoung]
    times = [told] + times + [tyoung]
    ranks = range(0, len(times))
    only_new = list()
    while missing > 0:
        if len(ranks) > 2:
            distrranks = list()
            for i in range(1, len(ranks)):
                temp = ranks[i] * (intp1(times[i-1], birth, death) - intp1(times[i], birth, death))
                distrranks.append(temp)
            try:
                distrranks = [x/sum(distrranks) for x in distrranks]
                for i in range(1, len(distrranks)):
                    distrranks[i] = distrranks[i] + distrranks[i-1]
                r = random.uniform(0, 1)
                addrank = min([idx for idx, x in enumerate(distrranks) if x > r])
            except ZeroDivisionError:
                addrank = 0
            except ValueError:
                addrank = 0
        else:
            addrank = 0
        r = random.uniform(0, 1)
        const = intp1(times[addrank], birth, death) - intp1(times[addrank+1], birth, death)
        try:
            temp = intp1(times[addrank+1], birth, death) / const
        except ZeroDivisionError:
            temp = 0.0
        xnew = 1 / (death - birth) * log((1 - (r + temp) * const * birth) / (1 - (r + temp) * const * death))
        ages.append(xnew)
        only_new.append(xnew)
        ages.sort(reverse=True)
        missing -= 1
    only_new.sort(reverse=True)
    return only_new


