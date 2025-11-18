"""Functions to handle various numerical operations, including optimization."""

import random
import sys
from collections.abc import Callable
from decimal import Decimal as D
from math import exp, log
from typing import Any

import numpy as np
from scipy.optimize import dual_annealing, minimize, minimize_scalar

# Raise on overflow
np.seterr(all="raise")


def get_bd(r: float, a: float) -> tuple[float, float]:
    """Convert turnover and relative extinction to birth and death rates.

    Args:
        r: Turnover rate (net diversification, birth - death).
        a: Relative extinction rate (death / birth).

    Returns:
        Tuple of (birth rate, death rate).
    """
    return -r / (a - 1), -a * r / (a - 1)


def get_ra(b: float, d: float) -> tuple[float, float]:
    """Convert birth and death rates to turnover and relative extinction.

    Args:
        b: Birth rate.
        d: Death rate.

    Returns:
        Tuple of (turnover rate, relative extinction rate).
    """
    return (b - d, d / b)


def wrapped_lik_constant(x: tuple[float, float], sampling: float, ages: list[float]) -> float:
    """Wrapper for birth-death likelihood function for optimization.

    Converts turnover and relative extinction parameters to birth/death rates
    before computing the likelihood.

    Args:
        x: Tuple of (turnover, relative extinction).
        sampling: Sampling fraction in (0, 1].
        ages: List of node ages.

    Returns:
        Negative log-likelihood value.
    """
    return lik_constant(get_bd(*x), sampling, ages)


def wrapped_lik_constant_yule(x: float, sampling: float, ages: list[float]) -> float:
    """Wrapper for Yule model likelihood function for optimization.

    Assumes zero extinction (pure birth process).

    Args:
        x: Birth rate.
        sampling: Sampling fraction in (0, 1].
        ages: List of node ages.

    Returns:
        Negative log-likelihood value.
    """
    return lik_constant((x, 0.0), sampling, ages)


def two_step_optim(
    func: Callable[..., float], x0: tuple[float, ...], bounds: tuple[tuple[float, float], ...], args: tuple[Any, ...]
) -> list[float]:
    """Optimize a function using a two-step approach.

    First attempts L-BFGS-B (fast gradient-based method), then falls back to
    simulated annealing if L-BFGS-B fails.

    Args:
        func: Objective function to minimize.
        x0: Initial parameter values.
        bounds: Parameter bounds as tuples of (min, max) for each parameter.
        args: Additional arguments to pass to the objective function.

    Returns:
        Optimized parameter values as a list.

    Raises:
        Exception: If both optimization methods fail.
    """
    try:
        result = minimize(func, x0=x0, bounds=bounds, args=args, method="L-BFGS-B")
        if result["success"]:
            return result["x"].tolist()  # type: ignore[no-any-return]
    except FloatingPointError:
        pass

    result = dual_annealing(func, x0=x0, bounds=bounds, args=args)
    if result["success"]:
        return result["x"].tolist()  # type: ignore[no-any-return]

    raise Exception(f"Optimization failed: {result['message']} (code {result['status']})")


def optim_bd(ages: list[float], sampling: float, min_bound: float = 1e-9) -> tuple[float, float]:
    """Optimize birth and death rates from node ages under a birth-death model.

    Uses maximum likelihood estimation with the Magallon-Sanderson crown estimator
    for initial values.

    Args:
        ages: List of node ages (splitting times).
        sampling: Sampling fraction in (0, 1].
        min_bound: Minimum allowed birth rate (default: 1e-9).

    Returns:
        Tuple of (optimized birth rate, optimized death rate).
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


def optim_yule(ages: list[float], sampling: float, min_bound: float = 1e-9) -> tuple[float, float]:
    """Optimize birth rate under a Yule (pure birth) model.

    Assumes zero extinction rate.

    Args:
        ages: List of node ages (splitting times).
        sampling: Sampling fraction in (0, 1].
        min_bound: Minimum allowed birth rate (default: 1e-9).

    Returns:
        Tuple of (optimized birth rate, 0.0).

    Raises:
        Exception: If optimization fails.
    """
    bounds = (min_bound, 100)
    result = minimize_scalar(wrapped_lik_constant_yule, bounds=bounds, args=(sampling, ages), method="Bounded")
    if result["success"]:
        return (result["x"], 0.0)

    raise Exception(f"Optimization failed: {result['message']} (code {result['status']})")


def p0_exact(t: float, l: float, m: float, rho: float) -> D:  # noqa: E741
    """Compute p0 using exact Decimal arithmetic.

    Used as fallback when floating-point arithmetic fails due to overflow.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).
        rho: Sampling fraction.

    Returns:
        Probability of no sampled descendants as a Decimal.
    """
    t_dec = D(t)
    l_dec = D(l)
    m_dec = D(m)
    rho_dec = D(rho)
    return D(1) - rho_dec * (l_dec - m_dec) / (
        rho_dec * l_dec + (l_dec * (D(1) - rho_dec) - m_dec) * (-(l_dec - m_dec) * t_dec).exp()
    )


def p0(t: float, l: float, m: float, rho: float) -> float:  # noqa: E741
    """Compute the probability of no sampled descendants.

    Probability that an individual alive at time `t` before present has no sampled
    descendants (extant or extinct), assuming no past sampling. Falls back to
    exact Decimal arithmetic if floating-point overflow occurs.

    Reference: Stadler (2010), Journal of Theoretical Biology 267(3):396-404,
    remark 3.2. Originally implemented as `TreePar:::p0`.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).
        rho: Sampling fraction.

    Returns:
        Probability of no sampled descendants.
    """
    try:
        return 1 - rho * (l - m) / (rho * l + (l * (1 - rho) - m) * exp(-(l - m) * t))
    except FloatingPointError:
        return float(p0_exact(t, l, m, rho))


def p1_exact(t: float, l: float, m: float, rho: float) -> D:  # noqa: E741
    """Compute p1 using exact Decimal arithmetic.

    Used as fallback when floating-point arithmetic fails due to overflow.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).
        rho: Sampling fraction.

    Returns:
        Probability of exactly one sampled descendant as a Decimal.
    """
    t_dec = D(t)
    l_dec = D(l)
    m_dec = D(m)
    rho_dec = D(rho)
    num = rho_dec * (l_dec - m_dec) ** D(2) * (-(l_dec - m_dec) * t_dec).exp()
    denom = (rho_dec * l_dec + (l_dec * (D(1) - rho_dec) - m_dec) * (-(l_dec - m_dec) * t_dec).exp()) ** D(2)
    return num / denom


def p1_orig(t: float, l: float, m: float, rho: float) -> float:  # noqa: E741
    """Original implementation of p1 for testing and comparison.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).
        rho: Sampling fraction.

    Returns:
        Probability of exactly one sampled descendant.
    """
    try:
        num = rho * (l - m) ** 2 * np.exp(-(l - m) * t)
        denom = (rho * l + (l * (1 - rho) - m) * np.exp(-(l - m) * t)) ** 2
        res: float = num / denom
    except (OverflowError, FloatingPointError):
        res = float(p1_exact(t, l, m, rho))
    if res == 0.0:
        return sys.float_info.min
    return res


def p1(t: float, l: float, m: float, rho: float) -> float:  # noqa: E741
    """Compute the probability of exactly one sampled descendant.

    Probability that an individual alive at time `t` before present has precisely
    one sampled extant descendant and no sampled extinct descendants, assuming
    no past sampling. Optimized version using common subexpression elimination.

    Reference: Stadler (2010), Journal of Theoretical Biology 267(3):396-404,
    remark 3.2. Originally implemented as `TreePar:::p1`.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).
        rho: Sampling fraction.

    Returns:
        Probability of exactly one sampled descendant.
    """
    try:
        ert = np.exp(-(l - m) * t, dtype=np.float64)
        num = rho * (l - m) ** 2 * ert
        denom = (rho * l + (l * (1 - rho) - m) * ert) ** 2
        res: float = num / denom
    except (OverflowError, FloatingPointError):
        res = float(p1_exact(t, l, m, rho))
    if res == 0.0:
        return sys.float_info.min
    return res


def intp1_exact(t: float, l: float, m: float) -> D:  # noqa: E741
    """Compute intp1 using exact Decimal arithmetic.

    Used as fallback when floating-point arithmetic fails due to overflow.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).

    Returns:
        Integration constant as a Decimal.
    """
    l_dec = D(l)
    m_dec = D(m)
    t_dec = D(t)
    num = D(1) - (-(l_dec - m_dec) * t_dec).exp()
    denom = l_dec - m_dec * (-(l_dec - m_dec) * t_dec).exp()
    return num / denom


def intp1(t: float, l: float, m: float) -> float:  # noqa: E741
    """Compute integration constant for sampling missing speciation event times.

    This is a portion of the cdf used to perform inverse-transform sampling of missing speciation event times
    under a constant-rate birth-death model. It is the c_2 term from equation A.2 in Cusimano et al. (2012),
    Systematic Biology 61(5):785-792. Originally implemented as `TreeSim:::intp1`.

    Args:
        t: Time before present.
        l: Birth rate (lambda).
        m: Death rate (mu).

    Returns:
        Integration constant value.
    """
    try:
        return (1 - exp(-(l - m) * t)) / (l - m * exp(-(l - m) * t))
    except OverflowError:
        return float(intp1_exact(t, l, m))


def lik_constant(
    vec: tuple[float, float],
    rho: float,
    t: list[float],
    root: int = 1,
    survival: int = 1,
    p1: Callable[[float, float, float, float], float] = p1,
) -> float:
    """Calculate likelihood of a constant-rate birth-death process.

    Likelihood conditioned on waiting times and incomplete sampling. Based on
    `TreePar::LikConstant` by Tanja Stadler. Reference: Stadler (2009), Journal
    of Theoretical Biology 261:58-66.

    Args:
        vec: Tuple of (birth rate, death rate).
        rho: Sampling fraction.
        t: List of waiting times (will be sorted in place).
        root: Include root contribution (1) or not (0). Default: 1.
        survival: Assume process survival (1) or not (0). Default: 1.
        p1: Function to compute p1 probability (default: `p1`).

    Returns:
        Negative log-likelihood value.
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


def crown_capture_probability(n: int, k: int) -> float:
    """Calculate probability of observing the crown node in an incomplete sample.

    Probability that a random sample of `k` taxa from a clade of `n` total taxa
    includes the crown (root) node, under a Yule process.

    Reference: Sanderson (1996), Systematic Biology 45:168-173.

    Args:
        n: Total number of taxa in the clade.
        k: Number of sampled taxa.

    Returns:
        Probability of crown node inclusion.

    Raises:
        Exception: If `n < k`.
    """
    if n < k:
        raise Exception(f"n must be greater than or equal to k (n={n}, k={k})")
    if n == 1 and k == 1:
        return 0.0  # not technically correct but it works for our purposes
    return 1 - 2 * (n - k) / ((n - 1) * (k + 1))


# TODO: This could probably be optimized
def get_new_times(
    ages: list[float],
    birth: float,
    death: float,
    missing: int,
    told: float | None = None,
    tyoung: float | None = None,
) -> list[float]:
    """Simulate new speciation event times in an incomplete phylogeny.

    Simulates missing speciation events under a constant-rate birth-death process.
    Adapted from `TreeSim::corsim` by Tanja Stadler. Reference: Cusimano et al.
    (2012), Systematic Biology 61(5):785-792.

    Args:
        ages: List of existing waiting times (will be sorted in place).
        birth: Birth rate.
        death: Death rate.
        missing: Number of missing taxa to simulate.
        told: Maximum simulated age. Defaults to `max(ages)`.
        tyoung: Minimum simulated age. Defaults to 0.

    Returns:
        List of simulated waiting times, sorted in descending order.

    Raises:
        Exception: If zero or negative branch lengths are detected.
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
    times = [told, *times, tyoung]
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
