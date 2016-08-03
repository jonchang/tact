# -*- coding: utf-8 -*-

from __future__ import division

import csv
import sys
import random
from math import log, exp, ceil
from decimal import Decimal as D

import dendropy
from scipy.optimize import minimize

def get_bd(r, a):
    r = float(r)
    a = float(a)
    b = r / (1 + a)
    d = b - r
    return b, d

def optim_bd(ages, sampling):
    """optimizes birth death using scipy"""
    return minimize(lambda x: lik_constant(x, sampling, ages), [1., 0.02], bounds=((sys.float_info.epsilon, None), (0, None)), method="TNC")["x"].tolist()


# from: TreePar::LikConstant
def lik_constant(vec, rho, t, root=1, survival=1):
    """vec: birth and death, rho: sampling, t: split times"""
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
    t = D(t)
    l = D(l)
    m = D(m)
    rho = D(rho)
    return rho*(l-m)**D(2) * (-(l-m)*t).exp()/(rho*l+(l*(1-rho)-m)*(-(l-m)*t).exp())**D(2)

def p1(t, l, m, rho):
    try:
        return rho*(l-m)**2 * exp(-(l-m)*t)/(rho*l+(l*(1-rho)-m)*exp(-(l-m)*t))**2
    except OverflowError:
        return float(p1_exact(t, l, m, rho))

def intp1_exact(t, l, m):
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

    This equation is taken from: Sanderson, M. J. 1996. How many taxa must
    be sampled to identify the root node of a large clade? Systematic Biology 45:168-173
    """
    return 1 - 2 * (n - k) / ((n - 1) * (k + 1))

def get_monophyletic_node(tree, species):
    mrca = tree.mrca(taxon_labels=species)
    if not mrca:
        return None
    new = set([x.taxon.label for x in mrca.leaf_iter()])
    if mrca and species.issuperset(new):
        return mrca

def get_birth_death_for_node(node, sampfrac):
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
    for x in node.preorder_internal_node_iter():
        if len(x.child_nodes()) != 2:
            return False
    return True


def get_short_branches(node):
    for edge in edge_iter(node):
        if edge.length <= 0.001:
            yield edge


def get_new_times(ages, birth, death, missing, told=None, tyoung=None):
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


