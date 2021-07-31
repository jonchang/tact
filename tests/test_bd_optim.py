from __future__ import division

from hypothesis import given, settings
import hypothesis.strategies as st

from tact.lib import optim_bd

# TODO: test with this array of ages for optimization edge cases:
#   st.just([38.428566000000004, 31.664333000000003, 29.987381, 27.317966000000002, 19.033624, 10.301403, 10.292682, 7.923767, 1.740758, 1.370161, 0.831707, 0.799889, 0.518386])

@settings(deadline=1500)
@given(st.lists(st.floats(min_value=0, allow_infinity=False, allow_nan=False, exclude_min=True), min_size=1), st.floats(min_value=1e-9, max_value=1))
def test_birth_death_can_run(ages, sampling):
    b, d = optim_bd(ages, sampling)
    assert b > 0
    assert d >= 0
    assert (b - d) > 0
    assert (d / b) <= 1


def test_birth_death(ages, sampling):
    optim_bd(ages, sampling)
