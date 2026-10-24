import hypothesis.strategies as st
from hypothesis import given

from tact.lib import crown_capture_probability


@st.composite
def clade_sizes(draw):
    n = draw(st.integers(min_value=2, max_value=10**6))
    k = draw(st.integers(min_value=1, max_value=n))
    return n, k


@given(clade_sizes())
def test_min_ccp_bounds_select_all_or_none(sizes):
    """--min-ccp 0 never applies, and --min-ccp 1 applies to every clade that is not fully sampled."""
    n, k = sizes
    ccp = crown_capture_probability(n, k)
    assert not ccp < 0
    assert (ccp < 1) == (k < n)
