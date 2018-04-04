from __future__ import division

from tact.lib import optim_bd_r, optim_bd_scipy, optim_bd_grid, optim_bd

from math import log, exp
import sys
import subprocess
import random

def get_new_times_r(ages, birth, death, missing, told=None):
    if told is None:
        told = max(ages)
    assert max(ages) <= told
    fmt = "source('corsim_new.R');cat(corsim_new(c({ages}), {birth}, {death}, {missing}, told={told}))".format(
            ages=",".join([str(x) for x in ages]),
            birth=birth,
            death=death,
            missing=missing,
            told=told
            )
    output = subprocess.check_output(["Rscript", "--vanilla", "--default-packages=base,stats", "-e", fmt])
    times = [float(x) for x in output.split()]
    return times

ages = [20.934955, 17.506532, 16.64467, 15.380987, 14.547092000000001, 13.664577999999999, 13.289480000000001, 11.667099, 9.799231, 9.510413, 9.029556000000001, 8.806255, 8.770727, 8.480102, 6.984475999999999, 6.706684, 2.11319, 0.545689, 0.147482]

sampling = 0.869565217391

def test_birth_death_scipy(benchmark):
    benchmark(optim_bd_scipy, ages, sampling)

#def test_birth_death_grid(benchmark):
#    benchmark(optim_bd_grid, ages, sampling)

def test_birth_death_mcmc(benchmark):
    benchmark(optim_bd, ages, sampling)

#def test_birth_death_r(benchmark):
#    benchmark(optim_bd_r, ages, sampling)


x = [21.115320000000004, 21.115320000000004, 18.234301000000002, 16.998473, 16.731151, 15.032467, 13.608031, 12.568441, 12.567588, 12.560298000000001, 12.070564000000001, 11.651428000000001, 10.873141, 10.66882, 10.459374, 9.880373, 9.71899, 9.686403, 9.527196, 9.036206, 8.247109, 7.696175, 7.428744, 6.204276, 6.189729, 5.463346, 4.846688, 4.54335, 4.096904, 2.533197, 2.505648, 2.381635, 0.46482599999999996, 0.375455, 0.321229, 0.00497]

birth = 0.174259100691
death = 0.0
missing = 78
told = max(x)
