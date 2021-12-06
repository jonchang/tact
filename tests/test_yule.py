from __future__ import division

from hypothesis import given, example
import hypothesis.strategies as st

from tact.lib import optim_yule


@example(
    ages=[
        5e-324,
        2.225073858507202e-308,
        0.3333333333333333,
        0.99999,
        1.1,
        746190.0,
        746190.0,
        746190.0,
        1.6377406994066752e16,
        1.7976931348623157e308,
    ],
    sampling=0.3477339723504501,
)
@given(
    st.lists(st.floats(min_value=0, allow_infinity=False, allow_nan=False, exclude_min=True), min_size=1),
    st.floats(min_value=1e-9, max_value=1),
)
def test_yule_has_no_extinction(ages, sampling):
    (b, d) = optim_yule(ages, sampling)
    assert d == 0
    assert b > 0
