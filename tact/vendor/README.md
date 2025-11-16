# Vendored Code

This directory contains vendored code from external projects, used to remove dependencies on scipy.

## Contents

- **pyprima/**: Pure Python implementation of COBYLA optimization algorithm
  - Source: scipy/_lib/pyprima/pyprima/
  - Git commit: 4899522c2a1a613b35bf5599680804ab641adb1d
  - License: BSD 3-Clause
  - See `pyprima/README.md` for details

- **scipy_optimize/**: Bounded scalar minimization function
  - Source: scipy/optimize/_optimize.py
  - Git commit: 4899522c2a1a613b35bf5599680804ab641adb1d
  - License: BSD 3-Clause (via SciPy)
  - See `scipy_optimize/README.md` for details

## Updating Vendored Code

Use the `update_vendor.sh` script in the repository root to update vendored code from a scipy repository:

```bash
./update_vendor.sh [path_to_scipy_repo]
```

If the path is not provided, it defaults to `../scipy`.

## License Compatibility

All vendored code uses BSD 3-Clause license, which is compatible with TACT's MIT license. When distributing TACT, you must:

1. Retain the original copyright notices
2. Include the BSD 3-Clause license text
3. Note that the code comes from SciPy/pyprima

See individual README files in each subdirectory for specific copyright information.

