"""Limit NumPy and SciPy to one linear algebra thread.

Import this before anything that imports NumPy.
"""

import os

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
