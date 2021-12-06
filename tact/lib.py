# -*- coding: utf-8 -*-

"""Functions to handle various numerical operations, including optimization."""

from __future__ import division

import random
import sys
from decimal import Decimal as D
from math import exp
from math import log

import numpy as np
from scipy.optimize import minimize, minimize_scalar, dual_annealing

# Raise on overflow
np.seterr(all="raise")


def get_bd(r, a):
    """
    Converts turnover and relative extinction to birth and death rates.

    Args:
        r (float): turnover or net diversification (birth - death)
        a (float): relative extinction (death / birth)

    Returns:
        (float, float): birth, death
    """
    return -r / (a - 1), -a * r / (a - 1)


def get_ra(b, d):
    """
    Converts birth and death to turnover and relative extinction rates.

    Args:
        b (float): birth rate
        d (float): extinction rate

    Returns:
        (float, float): turnover, relative extinction
    """
    return (b - d, d / b)


def wrapped_lik_constant(x, sampling, ages):
    """
    Wrapper for birth-death likelihood to make optimizing more convenient.

    Args:
        x (float, float): turnover, relative extinction
        sampling (float): sampling fraction (0, 1]
        ages (list): vector of node ages

    Returns:
        float: a likelihood
    """
    return lik_constant(get_bd(*x), sampling, ages)


def wrapped_lik_constant_yule(x, sampling, ages):
    """
    Wrapper for Yule likelihood to make optimizing more convenient.

    Args:
        x (float): birth rate
        sampling (float): sampling fraction (0, 1]
        ages (list): vector of node ages

    Returns:
        float: a likelihood
    """
    return lik_constant((x, 0.0), sampling, ages)


def two_step_optim(func, x0, bounds, args):
    """
    Conduct a two-step function optimization, first by using the fast L-BFGS-B method,
    and if that fails, use simulated annealing.

    Args:
        func (callable): function to optimize
        x0 (tuple): initial conditions
        bounds (tuple): boundary conditions
        args (lsit): additional argumnets to pass to `func`

    Returns:
        tuple: optimized parameter values
    """
    try:
        result = minimize(func, x0=x0, bounds=bounds, args=args, method="L-BFGS-B")
        if result["success"]:
            return result["x"].tolist()
    except FloatingPointError:
        pass

    result = dual_annealing(func, x0=x0, bounds=bounds, args=args)
    if result["success"]:
        return result["x"].tolist()

    raise Exception(f"Optimization failed: {result['message']} (code {result['status']})")


def optim_bd(ages, sampling, min_bound=1e-9):
    """
    Optimizes birth and death parameters given a vector of splitting times and sampling fraction.

    Args:
        ages (list): vector of node ages
        sampling (float): sampling fraction (0, 1]
        min_bound (float): minimum birth rate

    Returns:
        float, float: birth and death rates
    """
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
    """
    Optimizes birth parameter under a Yule model, given a vector of splitting times and sampling fraction.

    Args:
        ages (list): vector of node ages
        sampling (float): sampling fraction (0, 1]
        min_bound (float): minimum birth rate

    Returns:
        float, float: birth and death rates (where death is always 0)
    """
    bounds = (min_bound, 100)
    result = minimize_scalar(wrapped_lik_constant_yule, bounds=bounds, args=(sampling, ages), method="Bounded")
    if result["success"]:
        return (result["x"], 0.0)

    raise Exception(f"Optimization failed: {result['message']} (code {result['status']})")


def p0_exact(t, l, m, rho):  # noqa: E741
    "Exact version of `p0` using Decimal math."
    t = D(t)
    l = D(l)  # noqa: E741
    m = D(m)
    rho = D(rho)
    return D(1) - rho * (l - m) / (rho * l + (l * (D(1) - rho) - m) * (-(l - m) * t).exp())


def p0(t, l, m, rho):  # noqa: E741
    try:
        return 1 - rho * (l - m) / (rho * l + (l * (1 - rho) - m) * exp(-(l - m) * t))
    except FloatingPointError:
        return float(p0_exact(t, l, m, rho))


def p1_exact(t, l, m, rho):  # noqa: E741
    """Exact version of `p1` using Decimal math."""
    t = D(t)
    l = D(l)  # noqa: E741
    m = D(m)
    rho = D(rho)
    num = rho * (l - m) ** D(2) * (-(l - m) * t).exp()
    denom = (rho * l + (l * (1 - rho) - m) * (-(l - m) * t).exp()) ** D(2)
    return num / denom


def p1_orig(t, l, m, rho):  # noqa: E741
    """Original version of `p1`, here for testing and comparison purposes."""
    try:
        num = rho * (l - m) ** 2 * np.exp(-(l - m) * t)
        denom = (rho * l + (l * (1 - rho) - m) * np.exp(-(l - m) * t)) ** 2
        res = num / denom
    except (OverflowError, FloatingPointError):
        res = float(p1_exact(t, l, m, rho))
    if res == 0.0:
        return sys.float_info.min
    return res


def p1(t, l, m, rho):  # noqa: E741
    """
    Optimized version of `p1_orig` using common subexpression elimination and strength reduction
    from exponentiation to multiplication.
    """
    try:
        ert = np.exp(-(l - m) * t, dtype=np.float64)
        num = rho * (l - m) ** 2 * ert
        denom = (rho * l + (l * (1 - rho) - m) * ert) ** 2
        res = num / denom
    except (OverflowError, FloatingPointError):
        res = float(p1_exact(t, l, m, rho))
    if res == 0.0:
        return sys.float_info.min
    return res


def intp1_exact(t, l, m):  # noqa: E741
    """Exact version of `intp1` using Decimal math."""
    l = D(l)  # noqa: E741
    m = D(m)
    t = D(t)
    num = D(1) - (-(l - m) * t).exp()
    denom = l - m * (-(l - m) * t).exp()
    return num / denom


def intp1(t, l, m):  # noqa: E741
    try:
        return (1 - exp(-(l - m) * t)) / (l - m * exp(-(l - m) * t))
    except OverflowError:
        return float(intp1_exact(t, l, m))


def lik_constant(vec, rho, t, root=1, survival=1, p1=p1):
    """
    Calculates the likelihood of a constant-rate birth-death process, conditioned
    on the waiting times of a phylogenetic tree and degree of incomplete sampling.

    Based off of the R function `TreePar::LikConstant` written by Tanja Stadler.

    T. Stadler. On incomplete sampling under birth-death models and connections
    to the sampling-based coalescent. Jour. Theo. Biol. 261: 58-66, 2009.

    Args:
        vec (float, float): two element tuple of birth and death
        rho (float): sampling fraction
        t (list): vector of waiting times
        root (bool): include the root or not? (default: 1)
        survival (bool): assume survival of the process? (default: 1)

    Returns:
        float: a likelihood
    """
    l = vec[0]  # noqa: E741
    m = vec[1]
    t.sort(reverse=True)
    lik = (root + 1) * log(p1(t[0], l, m, rho))
    for tt in t[1:]:
        lik += log(l) + log(p1(tt, l, m, rho))
    if survival == 1:
        lik -= (root + 1) * log(1 - p0(t[0], l, m, rho))
    return -lik


def crown_capture_probability(n, k):
    """
    Calculate the probability that a sample of `k` taxa from a clade
    of `n` total taxa includes a root node, under a Yule process.

    This equation is taken from:

    Sanderson, M. J. 1996. How many taxa must be sampled to identify
    the root node of a large clade? Systematic Biology 45:168-173

    Args:
        n (int): total number of taxa
        k (int): sampled taxa

    Returns:
        float: probability
    """
    if n < k:
        raise Exception(f"n must be greater than or equal to k (n={n}, k={k})")
    if n == 1 and k == 1:
        return 0  # not technically correct but it works for our purposes
    return 1 - 2 * (n - k) / ((n - 1) * (k + 1))


# TODO: This could probably be optimized
def get_new_times(ages, birth, death, missing, told=None, tyoung=None):
    """
    Simulates new speciation events in an incomplete phylogeny assuming a
    constant-rate birth-death process.

    Adapted from the R function `TreeSim::corsim` written by Tanja Stadler.

    N. Cusimano, T. Stadler, S. Renner. A new method for handling missing
    species in diversification analysis applicable to randomly or
    non-randomly sampled phylogenies. Syst. Biol., 61(5): 785-792, 2012.

    Args:
        ages (list): vector of waiting times
        birth (float): birth rate
        death (float): death rate
        missing (int): number of missing taxa to simulate
        told (float): maximum simulated age (default: `max(ages)`)
        tyoung (float): minimum simulated age bound (default: `0`)

    Returns:
        list: vector of simulated waiting times.
    """
    if told is None:
        told = max(ages)
    if len(ages) > 0:
        if max(ages) > told and abs(max(ages) - told) > sys.float_info.epsilon:
            raise Exception("Zero or negative branch lengths detected in backbone phylogeny")
    if tyoung is None:
        tyoung = 0

    ages.sort(reverse=True)
    times = [x for x in ages if told >= x >= tyoung]
    times = [told] + times + [tyoung]
    ranks = range(0, len(times))
    only_new = []
    while missing > 0:
        if len(ranks) > 2:
            distrranks = []
            for i in range(1, len(ranks)):
                temp = ranks[i] * (intp1(times[i - 1], birth, death) - intp1(times[i], birth, death))
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
        xnew = 1 / (death - birth) * log((1 - (r + temp) * const * birth) / (1 - (r + temp) * const * death))
        only_new.append(xnew)
        missing -= 1
    only_new.sort(reverse=True)
    return only_new
