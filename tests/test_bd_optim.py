from __future__ import division

from hypothesis import given, settings
import hypothesis.strategies as st

from tact.lib import optim_bd


@settings(deadline=1500)
@given(st.lists(st.floats(min_value=0, allow_infinity=False, allow_nan=False, exclude_min=True), min_size=1), st.floats(min_value=1e-9, max_value=1))
def test_birth_death_can_run(x1, x2):
    b, d = optim_bd(x1, x2)
    assert b > 0
    assert d >= 0
    assert (b - d) > 0
    assert (d / b) <= 1


def test_birth_death(ages, sampling):
    optim_bd(ages, sampling)
