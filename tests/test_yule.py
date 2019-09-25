from __future__ import division

import pytest
from hypothesis import given, assume
import hypothesis.strategies as st

from tact.lib import optim_yule

@given(st.lists(st.floats(min_value=0, allow_infinity=False, allow_nan=False), min_size=1), st.floats(min_value=1e-9, max_value=1))
def test_yule_has_no_extinction(ages, sampling):
    (b, d) = optim_yule(ages, sampling)
    assert d == 0
    assert b > 0
