import hypothesis.strategies as st
import pytest
from hypothesis import assume, given

from tact.lib import lik_constant, p1, p1_exact, p1_orig


def test_lik_constant_exact(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1_exact)


def test_lik_constant_orig(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1_orig)


def test_lik_constant_opt(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1)


@given(
    birth=st.floats(min_value=1e-6, max_value=10),
    death=st.floats(min_value=0, max_value=10),
    sampling=st.floats(min_value=1e-9, max_value=1),
    ages=st.lists(st.floats(min_value=0, max_value=5000), min_size=1),
)
def test_lik_constants(birth, death, sampling, ages):
    assume(birth > death)
    assert lik_constant((birth, death), sampling, ages, p1=p1) == pytest.approx(
        lik_constant((birth, death), sampling, ages, p1=p1_orig)
    )
