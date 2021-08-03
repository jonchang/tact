from __future__ import division

import os
import pytest


@pytest.fixture
def datadir():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@pytest.fixture
def ages():
    return [20.934955, 17.506532, 16.64467, 15.380987, 14.547092000000001,
            13.664577999999999, 13.289480000000001, 11.667099, 9.799231,
            9.510413, 9.029556000000001, 8.806255, 8.770727, 8.480102,
            6.984475999999999, 6.706684, 2.11319, 0.545689, 0.147482]


@pytest.fixture
def sampling():
    return 0.869565217391


@pytest.fixture
def birth():
    return 0.09115578881894915


@pytest.fixture
def death():
    return 0.0
