# -*- coding: utf-8 -*-

# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

"""Functions ported from SciPy for optimization."""

from warnings import warn

from numpy import (atleast_1d, eye, argmin, zeros, shape, squeeze,
                   asarray, sqrt, Inf, asfarray)
import numpy as np


# standard status messages of optimizers
_status_message = {'success': 'Optimization terminated successfully.',
                   'maxfev': 'Maximum number of function evaluations has '
                              'been exceeded.',
                   'maxiter': 'Maximum number of iterations has been '
                              'exceeded.',
                   'pr_loss': 'Desired error not necessarily achieved due '
                              'to precision loss.',
                   'nan': 'NaN result encountered.',
                   'out_of_bounds': 'The result is outside of the provided '
                                    'bounds.'}


class OptimizeResult(dict):
    """ Represents the optimization result.
    Attributes
    ----------
    x : ndarray
        The solution of the optimization.
    success : bool
        Whether or not the optimizer exited successfully.
    status : int
        Termination status of the optimizer. Its value depends on the
        underlying solver. Refer to `message` for details.
    message : str
        Description of the cause of the termination.
    fun, jac, hess: ndarray
        Values of objective function, its Jacobian and its Hessian (if
        available). The Hessians may be approximations, see the documentation
        of the function in question.
    hess_inv : object
        Inverse of the objective function's Hessian; may be an approximation.
        Not available for all solvers. The type of this attribute may be
        either np.ndarray or scipy.sparse.linalg.LinearOperator.
    nfev, njev, nhev : int
        Number of evaluations of the objective functions and of its
        Jacobian and Hessian.
    nit : int
        Number of iterations performed by the optimizer.
    maxcv : float
        The maximum constraint violation.
    Notes
    -----
    `OptimizeResult` may have additional attributes not listed here depending
    on the specific solver being used. Since this class is essentially a
    subclass of dict with attribute accessors, one can see which
    attributes are available using the `OptimizeResult.keys` method.
    """

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


def is_array_scalar(x):
    """Test whether `x` is either a scalar or an array scalar.
    """
    return np.size(x) == 1


def _minimize_scalar_bounded(func, bounds, args=(),
                             xatol=1e-5, maxiter=500,
                             **unknown_options):
    """
    Options
    -------
    maxiter : int
        Maximum number of iterations to perform.
    xatol : float
        Absolute error in solution `xopt` acceptable for convergence.
    """
    maxfun = maxiter
    # Test bounds are of correct form
    if len(bounds) != 2:
        raise ValueError('bounds must have two elements.')
    x1, x2 = bounds

    if not (is_array_scalar(x1) and is_array_scalar(x2)):
        raise ValueError("Optimization bounds must be scalars"
                         " or array scalars.")
    if x1 > x2:
        raise ValueError("The lower bound exceeds the upper bound.")

    flag = 0

    sqrt_eps = sqrt(2.2e-16)
    golden_mean = 0.5 * (3.0 - sqrt(5.0))
    a, b = x1, x2
    fulc = a + golden_mean * (b - a)
    nfc, xf = fulc, fulc
    rat = e = 0.0
    x = xf
    fx = func(x, *args)
    num = 1
    fu = np.inf

    ffulc = fnfc = fx
    xm = 0.5 * (a + b)
    tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0
    tol2 = 2.0 * tol1

    while (np.abs(xf - xm) > (tol2 - 0.5 * (b - a))):
        golden = 1
        # Check for parabolic fit
        if np.abs(e) > tol1:
            golden = 0
            r = (xf - nfc) * (fx - ffulc)
            q = (xf - fulc) * (fx - fnfc)
            p = (xf - fulc) * q - (xf - nfc) * r
            q = 2.0 * (q - r)
            if q > 0.0:
                p = -p
            q = np.abs(q)
            r = e
            e = rat

            # Check for acceptability of parabola
            if ((np.abs(p) < np.abs(0.5*q*r)) and (p > q*(a - xf)) and
                    (p < q * (b - xf))):
                rat = (p + 0.0) / q
                x = xf + rat

                if ((x - a) < tol2) or ((b - x) < tol2):
                    si = np.sign(xm - xf) + ((xm - xf) == 0)
                    rat = tol1 * si
            else:      # do a golden-section step
                golden = 1

        if golden:  # do a golden-section step
            if xf >= xm:
                e = a - xf
            else:
                e = b - xf
            rat = golden_mean*e

        si = np.sign(rat) + (rat == 0)
        x = xf + si * np.maximum(np.abs(rat), tol1)
        fu = func(x, *args)
        num += 1

        if fu <= fx:
            if x >= xf:
                a = xf
            else:
                b = xf
            fulc, ffulc = nfc, fnfc
            nfc, fnfc = xf, fx
            xf, fx = x, fu
        else:
            if x < xf:
                a = x
            else:
                b = x
            if (fu <= fnfc) or (nfc == xf):
                fulc, ffulc = nfc, fnfc
                nfc, fnfc = x, fu
            elif (fu <= ffulc) or (fulc == xf) or (fulc == nfc):
                fulc, ffulc = x, fu

        xm = 0.5 * (a + b)
        tol1 = sqrt_eps * np.abs(xf) + xatol / 3.0
        tol2 = 2.0 * tol1

        if num >= maxfun:
            flag = 1
            break

    if np.isnan(xf) or np.isnan(fx) or np.isnan(fu):
        flag = 2

    fval = fx

    result = OptimizeResult(fun=fval, status=flag, success=(flag == 0),
                            message={0: 'Solution found.',
                                     1: 'Maximum number of function calls '
                                        'reached.',
                                     2: _status_message['nan']}.get(flag, ''),
                            x=xf, nfev=num, nit=num)

    return result


def minimize_scalar(fun, bracket=None, bounds=None, args=(),
                    method='brent', tol=None, options=None):
    """Minimization of scalar function of one variable.
    Parameters
    ----------
    fun : callable
        Objective function.
        Scalar function, must return a scalar.
    bracket : sequence, optional
        For methods 'brent' and 'golden', `bracket` defines the bracketing
        interval and can either have three items ``(a, b, c)`` so that
        ``a < b < c`` and ``fun(b) < fun(a), fun(c)`` or two items ``a`` and
        ``c`` which are assumed to be a starting interval for a downhill
        bracket search (see `bracket`); it doesn't always mean that the
        obtained solution will satisfy ``a <= x <= c``.
    bounds : sequence, optional
        For method 'bounded', `bounds` is mandatory and must have two items
        corresponding to the optimization bounds.
    args : tuple, optional
        Extra arguments passed to the objective function.
    method : str or callable, optional
        Type of solver.  Should be one of:
            - :ref:`Brent <optimize.minimize_scalar-brent>`
            - :ref:`Bounded <optimize.minimize_scalar-bounded>`
            - :ref:`Golden <optimize.minimize_scalar-golden>`
            - custom - a callable object (added in version 0.14.0), see below
        See the 'Notes' section for details of each solver.
    tol : float, optional
        Tolerance for termination. For detailed control, use solver-specific
        options.
    options : dict, optional
        A dictionary of solver options.
            maxiter : int
                Maximum number of iterations to perform.
            disp : bool
                Set to True to print convergence messages.
        See :func:`show_options()` for solver-specific options.
    Returns
    -------
    res : OptimizeResult
        The optimization result represented as a ``OptimizeResult`` object.
        Important attributes are: ``x`` the solution array, ``success`` a
        Boolean flag indicating if the optimizer exited successfully and
        ``message`` which describes the cause of the termination. See
        `OptimizeResult` for a description of other attributes.
    See also
    --------
    minimize : Interface to minimization algorithms for scalar multivariate
        functions
    show_options : Additional options accepted by the solvers
    Notes
    -----
    This section describes the available solvers that can be selected by the
    'method' parameter. The default method is *Brent*.
    Method :ref:`Brent <optimize.minimize_scalar-brent>` uses Brent's
    algorithm to find a local minimum.  The algorithm uses inverse
    parabolic interpolation when possible to speed up convergence of
    the golden section method.
    Method :ref:`Golden <optimize.minimize_scalar-golden>` uses the
    golden section search technique. It uses analog of the bisection
    method to decrease the bracketed interval. It is usually
    preferable to use the *Brent* method.
    Method :ref:`Bounded <optimize.minimize_scalar-bounded>` can
    perform bounded minimization. It uses the Brent method to find a
    local minimum in the interval x1 < xopt < x2.
    **Custom minimizers**
    It may be useful to pass a custom minimization method, for example
    when using some library frontend to minimize_scalar. You can simply
    pass a callable as the ``method`` parameter.
    The callable is called as ``method(fun, args, **kwargs, **options)``
    where ``kwargs`` corresponds to any other parameters passed to `minimize`
    (such as `bracket`, `tol`, etc.), except the `options` dict, which has
    its contents also passed as `method` parameters pair by pair.  The method
    shall return an `OptimizeResult` object.
    The provided `method` callable must be able to accept (and possibly ignore)
    arbitrary parameters; the set of parameters accepted by `minimize` may
    expand in future versions and then these parameters will be passed to
    the method. You can find an example in the scipy.optimize tutorial.
    .. versionadded:: 0.11.0
    Examples
    --------
    Consider the problem of minimizing the following function.
    >>> def f(x):
    ...     return (x - 2) * x * (x + 2)**2
    Using the *Brent* method, we find the local minimum as:
    >>> from scipy.optimize import minimize_scalar
    >>> res = minimize_scalar(f)
    >>> res.x
    1.28077640403
    Using the *Bounded* method, we find a local minimum with specified
    bounds as:
    >>> res = minimize_scalar(f, bounds=(-3, -1), method='bounded')
    >>> res.x
    -2.0000002026
    """
    if not isinstance(args, tuple):
        args = (args,)

    if callable(method):
        meth = "_custom"
    else:
        meth = method.lower()
    if options is None:
        options = {}

    if tol is not None:
        options = dict(options)
        if meth == 'bounded' and 'xatol' not in options:
            warn("Method 'bounded' does not support relative tolerance in x; "
                 "defaulting to absolute tolerance.", RuntimeWarning)
            options['xatol'] = tol
        elif meth == '_custom':
            options.setdefault('tol', tol)
        else:
            options.setdefault('xtol', tol)

    # replace boolean "disp" option, if specified, by an integer value.
    disp = options.get('disp')
    if isinstance(disp, bool):
        options['disp'] = 2 * int(disp)

    if meth == '_custom':
        return method(fun, args=args, bracket=bracket, bounds=bounds, **options)
    elif meth == 'bounded':
        if bounds is None:
            raise ValueError('The `bounds` parameter is mandatory for '
                             'method `bounded`.')
        return _minimize_scalar_bounded(fun, bounds, args, **options)
    else:
        raise ValueError('Unknown solver %s' % method)
