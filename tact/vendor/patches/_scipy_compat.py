"""
Compatibility layer for scipy.optimize classes used by pyprima.

This module provides minimal implementations of scipy.optimize classes
needed by pyprima, allowing it to work without scipy as a dependency.
"""

import numpy as np


class Bounds:
    """Minimal Bounds class compatible with scipy.optimize.Bounds."""
    
    def __init__(self, lb, ub, keep_feasible=False):
        """
        Parameters
        ----------
        lb, ub : array_like, optional
            Lower and upper bounds on independent variables. Each array must
            have the same size as x or be a scalar, in which case a bound will
            be the same for all variables. Use np.inf with an appropriate sign
            to disable bounds on all or some variables.
        keep_feasible : bool, optional
            Whether to keep the solution within the bounds.
        """
        self.lb = np.asarray(lb, dtype=np.float64)
        self.ub = np.asarray(ub, dtype=np.float64)
        self.keep_feasible = keep_feasible


class LinearConstraint:
    """Minimal LinearConstraint class compatible with scipy.optimize.LinearConstraint."""
    
    def __init__(self, A, lb, ub, keep_feasible=False):
        """
        Parameters
        ----------
        A : {array_like, sparse matrix}, shape (m, n)
            Matrix defining the constraint.
        lb, ub : array_like
            Lower and upper bounds on the constraint. Each array must have the
            same size as the number of rows of A or be a scalar, in which case
            a bound will be the same for all constraints. Use np.inf with an
            appropriate sign to disable bounds on all or some constraints.
        keep_feasible : bool, optional
            Whether to keep the solution within the bounds.
        """
        self.A = np.asarray(A)
        self.lb = np.asarray(lb, dtype=np.float64)
        self.ub = np.asarray(ub, dtype=np.float64)
        self.keep_feasible = keep_feasible


class NonlinearConstraint:
    """Minimal NonlinearConstraint class compatible with scipy.optimize.NonlinearConstraint."""
    
    def __init__(self, fun, lb, ub, jac=None, hess=None, keep_feasible=False):
        """
        Parameters
        ----------
        fun : callable
            The function defining the constraint.
        lb, ub : array_like
            Lower and upper bounds on the constraint. Each array must have the
            same size as the number of constraints or be a scalar, in which case
            a bound will be the same for all constraints. Use np.inf with an
            appropriate sign to disable bounds on all or some constraints.
        jac : {callable, '2-point', '3-point', 'cs'}, optional
            Method of computing the Jacobian matrix.
        hess : {callable, '2-point', '3-point', 'cs', HessianUpdateStrategy}, optional
            Method of computing the Hessian matrix.
        keep_feasible : bool, optional
            Whether to keep the solution within the bounds.
        """
        self.fun = fun
        self.lb = np.asarray(lb, dtype=np.float64)
        self.ub = np.asarray(ub, dtype=np.float64)
        self.jac = jac
        self.hess = hess
        self.keep_feasible = keep_feasible


class OptimizeResult(dict):
    """Minimal OptimizeResult class compatible with scipy.optimize.OptimizeResult."""
    
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            return f"{self.__class__.__name__}({dict.__repr__(self)})"
        else:
            return self.__class__.__name__ + "()"

