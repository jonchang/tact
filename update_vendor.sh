#!/bin/bash
# Script to update vendored code from scipy repository
#
# Usage: ./update_vendor.sh [scipy_repo_path]
#
# If scipy_repo_path is not provided, defaults to ../scipy

set -e

SCIPY_REPO="${1:-../scipy}"

if [ ! -d "$SCIPY_REPO" ]; then
    echo "Error: SciPy repository not found at $SCIPY_REPO"
    echo "Usage: $0 [scipy_repo_path]"
    exit 1
fi

# Get the git commit hash
COMMIT=$(cd "$SCIPY_REPO" && git rev-parse HEAD)
echo ">>> Updating vendored code from SciPy commit: $COMMIT"

# Update pyprima
echo ">>> Updating pyprima..."
rm -rf tact/vendor/pyprima/*
for item in "$SCIPY_REPO/scipy/_lib/pyprima/pyprima/"*; do
    cp -r "$item" tact/vendor/pyprima/
done
# Copy license file separately
cp "$SCIPY_REPO/scipy/_lib/pyprima/LICENCE.txt" tact/vendor/pyprima/

# Drop bundled examples which are not used at runtime
rm -rf tact/vendor/pyprima/src/pyprima/examples
rm -rf tact/vendor/pyprima/tests
rm -rf tact/vendor/pyprima/pyproject.toml

# Write our custom README file
cat <<EOF > tact/vendor/pyprima/README.md
# Vendored pyprima code

This directory contains vendored code from pyprima, specifically the \`cobyla\` solver.

## Source

- **Repository**: https://github.com/scipy/scipy
- **Path**: \`scipy/_lib/pyprima/pyprima/\`
- **Git commit**: $COMMIT
- **License**: BSD 3-Clause (see [LICENSE.txt](LICENSE.txt))

## About

Refer to the original README file for more information: https://github.com/scipy/scipy/blob/main/scipy/_lib/pyprima/pyprima/README.md

## Updating

Use the \`update_vendor.sh\` script in the repository root to update this vendored code.
EOF

# Apply local patches to restore scipy compatibility shim
echo ""
echo ">>> Applying local patches..."
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PATCHES_DIR="$SCRIPT_DIR/tact/vendor/patches"
PYPRIMA_DIR="$SCRIPT_DIR/tact/vendor/pyprima"

# Copy the compatibility shim file
if [ -f "$PATCHES_DIR/_scipy_compat.py" ]; then
    cp "$PATCHES_DIR/_scipy_compat.py" "$PYPRIMA_DIR/src/pyprima/_scipy_compat.py"
    echo "Copied _scipy_compat.py"
else
    echo "Warning: _scipy_compat.py not found in $PATCHES_DIR"
fi

# Apply patches
cd "$PYPRIMA_DIR"
for patch_file in "$PATCHES_DIR"/*.patch; do
    patch -p1 < "$patch_file"
done
cd - > /dev/null

echo ""
echo ">>> Update complete!"

echo "Note: The _minimize_scalar_bounded function is manually extracted from scipy/optimize/_optimize.py"
echo "      and placed in tact/vendor/scipy_optimize/_minimize_scalar_bounded.py."
