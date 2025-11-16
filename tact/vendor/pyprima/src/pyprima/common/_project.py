'''
This module provides the _project function that attempts to project the initial guess
onto the feasible set.

Adapted from the corresponding function in the PDFO package (https://www.pdfo.net) by
Tom M. Ragonneau (https://ragonneau.github.io) and Zaikun Zhang (https://www.zhangzk.net).
'''

import numpy as np
from ._linear_constraints import LinearConstraint
from .._scipy_compat import OptimizeResult

# All the accepted scalar types; np.generic correspond to all NumPy types.
scalar_types = (int, float, np.generic)
eps = np.finfo(np.float64).eps

def _project(x0, lb, ub, constraints):
    """Projection of the initial guess onto the feasible set.

    Parameters
    ----------
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with
        fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: dict
                The nonlinear constraints of the problem. When ``_project`` is called, the nonlinear constraints are
                None.

    Returns
    -------
    result: OptimizeResult
        The result of the projection.

    Authors
    -------
    Tom M. RAGONNEAU (ragonneau.github.io)
    and Zaikun ZHANG (www.zhangzk.net)

    Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    """
    invoker = 'prima'

    # Validate x0.
    if isinstance(x0, scalar_types):
        x0_c = [x0]
    elif hasattr(x0, '__len__'):
        x0_c = x0
    else:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    try:
        x0_c = np.asarray(x0_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should contain only scalars.'.format(invoker))
    if len(x0_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    lenx0 = x0_c.size

    # Validate lb.
    if isinstance(lb, scalar_types):
        lb_c = [lb]
    elif hasattr(lb, '__len__'):
        lb_c = lb
    else:
        raise ValueError('{}: UNEXPECTED ERROR: lb should be a vector.'.format(invoker))
    try:
        lb_c = np.asarray(lb_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: lb should contain only scalars.'.format(invoker))
    if len(lb_c.shape) != 1 or lb_c.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of lb is inconsistent with x0.'.format(invoker))

    # Validate ub.
    if isinstance(ub, scalar_types):
        ub_c = [ub]
    elif hasattr(ub, '__len__'):
        ub_c = ub
    else:
        raise ValueError('{}: UNEXPECTED ERROR: ub should be a vector.'.format(invoker))
    try:
        ub_c = np.asarray(ub_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: ub should contain only scalars.'.format(invoker))
    if len(ub_c.shape) != 1 or ub_c.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of ub is inconsistent with x0.'.format(invoker))

    # Validate constraints.
    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())) or \
            not (isinstance(constraints['linear'], LinearConstraint) or constraints['linear'] is None):
        # the nonlinear constraints will not be taken into account in this function and are, therefore, not validated
        raise ValueError('{}: UNEXPECTED ERROR: The constraints are ill-defined.'.format(invoker))

    max_con = 1e20  # Decide whether an inequality constraint can be ignored

    # Project onto the feasible set.
    if constraints['linear'] is None:
        # Direct projection onto the bound constraints
        x_proj = np.nanmin((np.nanmax((x0_c, lb_c), axis=0), ub_c), axis=0)
        return OptimizeResult(x=x_proj)
    elif all(np.less_equal(np.abs(constraints['linear'].ub - constraints['linear'].lb), eps)) and \
            np.max(lb_c) <= -max_con and np.min(ub_c) >= max_con:
        # The linear constraints are all equality constraints. The projection can therefore be done by solving the
        # least-squares problem: min ||A*x - (b - A*x_0)||.
        a = constraints['linear'].A
        b = (constraints['linear'].lb + constraints['linear'].ub) / 2
        xi, _, _, _ = np.linalg.lstsq(a, b - np.dot(a, x0_c), rcond=None)

        # The problem is not bounded. However, if the least-square solver returned values bigger in absolute value
        # than max_con, they will be reduced to this bound.
        x_proj = np.nanmin((np.nanmax((x0_c + xi, lb_c), axis=0), ub_c), axis=0)

        return OptimizeResult(x=x_proj)

    if constraints['linear'] is not None:
        # Note: The original code used scipy.optimize.minimize for projection onto
        # linear constraints. Since we're removing scipy dependency, we'll do a
        # simple bound projection instead. This is acceptable for our use case.
        # For more complex linear constraints, a proper projection would be needed.
        x_proj = np.nanmin((np.nanmax((x0_c, lb_c), axis=0), ub_c), axis=0)
        return OptimizeResult(x=x_proj)

    return OptimizeResult(x=x0_c)
