"""
Vendored code from scipy.optimize._optimize

This module contains the _minimize_scalar_bounded function from SciPy,
vendored to remove the scipy dependency.

Source: scipy/optimize/_optimize.py
Git commit: 4899522c2a1a613b35bf5599680804ab641adb1d
License: BSD 3-Clause (via SciPy)
Copyright: SciPy Developers
"""

from ._minimize_scalar_bounded import minimize_scalar_bounded

__all__ = ['minimize_scalar_bounded']

