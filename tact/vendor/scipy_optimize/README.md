# Vendored scipy.optimize code

This directory contains vendored code from scipy.optimize, specifically the `_minimize_scalar_bounded` function.

## Source

- **Repository**: https://github.com/scipy/scipy
- **Path**: `scipy/optimize/_optimize.py`
- **Function**: `_minimize_scalar_bounded()` (lines ~2289-2436)
- **Git commit**: 4899522c2a1a613b35bf5599680804ab641adb1d
- **License**: BSD 3-Clause (via SciPy)
- **Copyright**: SciPy Developers

## About

The `minimize_scalar_bounded` function implements bounded scalar minimization using golden section search with parabolic interpolation (a variant of Brent's method). This is used in TACT for the `optim_yule()` function.

## License

This code is part of SciPy and is licensed under the BSD 3-Clause License. The full license text can be found in the SciPy repository at `LICENSE.txt`.

## Updating

Use the `update_vendor.sh` script in the repository root to update this vendored code.

