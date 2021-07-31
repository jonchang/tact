from __future__ import division

from hypothesis import given, settings, example
import hypothesis.strategies as st

from tact.lib import optim_bd

@settings(deadline=2000)
@example(ages=[38.428566000000004, 31.664333000000003, 29.987381, 27.317966000000002, 19.033624, 10.301403, 10.292682, 7.923767, 1.740758, 1.370161, 0.831707, 0.799889, 0.518386], sampling=0.9931334704375381)
@example(ages=[5e-324, 2.00001, 1.7976931348623157e+308, 0.5, 655944740946.0, 626147432795.0, 0.3333333333333333, 8.881784197001251e-16, 8.881784197001251e-16], sampling=0.06298810337307759)
@given(st.lists(st.floats(min_value=0, allow_infinity=False, allow_nan=False, exclude_min=True), min_size=1), st.floats(min_value=1e-9, max_value=1))
def test_birth_death_can_run(ages, sampling):
    b, d = optim_bd(ages, sampling)
    assert b > 0
    assert d >= 0
    assert (b - d) > 0
    assert (d / b) <= 1


def test_birth_death(ages, sampling):
    optim_bd(ages, sampling)
