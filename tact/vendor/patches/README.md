# Vendor Patches

This directory contains patches and compatibility shims for the vendored pyprima code.

## Contents

- **`_scipy_compat.py`**: Compatibility shim that provides minimal implementations of scipy.optimize classes (`Bounds`, `LinearConstraint`, `NonlinearConstraint`, `OptimizeResult`) to allow pyprima to work without scipy as a dependency.

- **`*.patch`**: Unified diff patches that modify the vendored pyprima code to use the compatibility shim instead of importing from scipy directly.

## Patches

1. **`01-__init__.py.patch`**: Changes imports in `src/pyprima/__init__.py` to use the compatibility shim.
2. **`02-_bounds.py.patch`**: Changes imports in `src/pyprima/common/_bounds.py` to use the compatibility shim.
3. **`03-_linear_constraints.py.patch`**: Changes imports in `src/pyprima/common/_linear_constraints.py` to use the compatibility shim.
4. **`04-_project.py.patch`**: Changes imports in `src/pyprima/common/_project.py` to use the compatibility shim.

## Application

These patches are automatically applied by the `update_vendor.sh` script after updating the vendored code from the scipy repository. The script:

1. Copies `_scipy_compat.py` to the pyprima source directory
2. Applies all `.patch` files in order

## Manual Application

If you need to apply patches manually:

```bash
cd tact/vendor/pyprima
cp ../patches/_scipy_compat.py src/pyprima/
patch -p1 < ../patches/01-__init__.py.patch
patch -p1 < ../patches/02-_bounds.py.patch
patch -p1 < ../patches/03-_linear_constraints.py.patch
patch -p1 < ../patches/04-_project.py.patch
```

