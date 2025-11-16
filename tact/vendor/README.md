# Vendored Code

This directory contains vendored code from external projects, used to remove dependencies on scipy.

## Contents

- [**pyprima/**](pyprima/README.md): Pure Python implementation of COBYLA optimization algorithm
- [**scipy_optimize/**](scipy_optimize/README.md): Bounded scalar minimization function
- [**patches/**](patches/README.md): Various patches needed for the vendored code to work

## Updating Vendored Code

Use the `update_vendor.sh` script in the repository root to update vendored code from a scipy repository.

## License Compatibility

All vendored code uses BSD 3-Clause license, which is compatible with TACT's MIT license. When distributing TACT, you must:

1. Retain the original copyright notices
2. Include the BSD 3-Clause license text
3. Note that the code comes from SciPy/pyprima

See individual README files in each subdirectory for specific copyright information.
