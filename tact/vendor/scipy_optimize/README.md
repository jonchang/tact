# Vendored scipy.optimize code

This directory contains vendored code from scipy.optimize, specifically the `_minimize_scalar_bounded` function.

## Source

- **Repository**: https://github.com/scipy/scipy
- **Path**: `scipy/optimize/_optimize.py`
- **Function**: `_minimize_scalar_bounded()` (lines ~2289-2436)
- **Git commit**: 4899522c2a1a613b35bf5599680804ab641adb1d
- **License**: BSD 3-Clause (see [`LICENSE.txt`](LICENSE.txt))

## About

The `minimize_scalar_bounded` function implements bounded scalar minimization using golden section search with parabolic interpolation (a variant of Brent's method). This is used in TACT for the `optim_yule()` function.
