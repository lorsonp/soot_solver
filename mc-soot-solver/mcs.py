"""This is a Stochastic Monte Carlo solver for the particle size distribution of soot."""

import math
import numpy as np
import random
import matplotlib.pyplot as plt

# ======================== #
# Generate particle system #
# ======================== #
"""Import source terms
Import initial concentrations of C16H10, C2H2, H, H2, O2, and OH.
If applicable, import initial particle distribution"""
n = []
start_time = 0.0
stop_time = 10.0


def sample_times(start_time, stop_time):
    sample_times = np.arange(start_time, stop_time, 0.05)
    return sample_times

# =============== #
# Calculate rates #
# =============== #

def inception_rate(reaction_temp, m_pyrene, pyrene_concentration):
    dA = 1.395*(3.**0.5)
    Boltzmann_k = 1.38064852*10**-23  #m^2 kg s^-2 K^-1
    pyrene_coagulation = (17.6/3.)*(dA**2)*m_pyrene*(np.pi*Boltzmann_k*reaction_temp/(2*m_pyrene))**0.5
    rate_of_inception = 0.5*pyrene_coagulation*pyrene_concentration**2.
    return rate_of_inception


def coagulation_rate(N, n):
    kernel_1 = 0
    kernel_2a = 0
    kernel_2b = 0
    for index in range(0, len(n)-1):
        kernel_1 += n(index)**(1./6.)
        kernel_2a += n(index)**(2./3.)
        kernel_2b += n(index)**(-1./2.)
    rate_of_coagulation = (1.4178/N)*((len(n)-2)*kernel_1 + kernel_2a*kernel_2b)
    return rate_of_coagulation


# =========================================== #
# Wait an exponentially distributed time step #
# =========================================== #


def time_step(inception_rate, coagulation_rate):
    r = random.uniform(0, 1)
    waiting_parameter = inception_rate + coagulation_rate
    tao = waiting_parameter**(-1)*math.log(1/r)
    return tao


# ============================== #
# Select event probabilistically #
# ============================== #


def select_event(inception_rate, coagulation_rate):
    r = random.uniform(0, 1)
    if r >= 0 and r <= (inception_rate/(inception_rate + coagulation_rate)):
        selection = 1
    else:
        selection = 2
    return selection


# ============== #
# Inception Step #
# ============== #


def inception_step(n):
    n = n.append(32)
    return n


"""def coagulation_step()
"""

test = inception_rate(1900, .0000000000000000000001, 16)
test2 = coagulation_rate(100, n)
test3 = select_event(test, test2)
print(test3)