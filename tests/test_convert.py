import hypothesis.strategies as st
import pytest
from hypothesis import assume, given

from tact.lib import get_bd, get_ra


@given(
    st.floats(min_value=1e-6, max_value=1e3, allow_infinity=False, allow_nan=False),
    st.floats(min_value=0, max_value=1 - 1e-6, allow_infinity=False, allow_nan=False),
)
def test_ra_inverts_bd(r, a):
    assert (r, a) == pytest.approx(get_ra(*get_bd(r, a)))


@given(
    st.floats(min_value=1e-6, max_value=1e3, allow_infinity=False, allow_nan=False),
    st.floats(min_value=0, allow_infinity=False, allow_nan=False),
)
def test_bd_inverts_ra(b, d):
    assume(b > d)
    assert (b, d) == pytest.approx(get_bd(*get_ra(b, d)))
