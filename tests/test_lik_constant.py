from __future__ import division
import pytest

from tact.lib import p1, p1_orig, p1_exact, lik_constant

def test_lik_constant_exact(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1_exact)

def test_lik_constant_orig(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1_orig)

def test_lik_constant_opt(benchmark, birth, death, sampling, ages):
    benchmark(lik_constant, (birth, death), sampling, ages, p1=p1)

def test_lik_constants(birth, death, sampling, ages):
    assert lik_constant((birth, death), sampling, ages, p1=p1) == pytest.approx(lik_constant((birth, death), sampling, ages, p1=p1_orig))
