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

x = [21.115320000000004, 21.115320000000004, 18.234301000000002, 16.998473, 16.731151, 15.032467, 13.608031, 12.568441, 12.567588, 12.560298000000001, 12.070564000000001, 11.651428000000001, 10.873141, 10.66882, 10.459374, 9.880373, 9.71899, 9.686403, 9.527196, 9.036206, 8.247109, 7.696175, 7.428744, 6.204276, 6.189729, 5.463346, 4.846688, 4.54335, 4.096904, 2.533197, 2.505648, 2.381635, 0.46482599999999996, 0.375455, 0.321229, 0.00497]

missing = 78
told = max(x)
