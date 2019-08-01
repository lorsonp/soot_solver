"""This is a Stochastic Monte Carlo solver for the particle size distribution of soot."""

import math
import numpy as np
import random
import cantera as ct
import matplotlib.pyplot as plt

# ======================== #
# Generate particle system #
# ======================== #
"""Import source terms:
- Need Temperature and pyrene concentration.
- If applicable, import initial particle distribution"""

gas1 = ct.Solution('gri30.xml')

# initiate matrix to track particles
# - each row is a separate 'group'
# - column index + 32 indicates the size of the particle (# of carbon atoms)
# - index values denote the number of particles of that size
# particles = np.zeros((100, 100))

N = 100  # sample size

start_time = 0.0
stop_time = 10.0


def sample_times(start_time, stop_time):
    sample_times = np.arange(start_time, stop_time, 0.05)
    return sample_times

# =============== #
# Calculate rates #
# =============== #

def inception_rate(reaction_temp, pyrene_concentration):
    """
    This function calculates the inception rate.
        Inputs:
            reaction_temp = temperature of the flame
            pyrene_concentration = # of pyrene molecules/total volume
        Outputs:
            The rate of inception (1/(m^3*s))
    """
    dA = 1.395e-10*(3.**0.5) # m
    reduced_m_pyrene = (202.256e-3/(6.0221409e23))*(0.5) # kg
    Boltzmann_k = 1.38064852e-23  #m^2 kg s^-2 K^-1
    pyrene_coagulation = (2.2*2/3)*dA**2*(16)*(np.pi*Boltzmann_k*reaction_temp/(2*reduced_m_pyrene))**0.5
    rate_of_inception = 0.5*pyrene_coagulation*pyrene_concentration**2.
    return rate_of_inception


def coagulation_rate(N, particles):
    """
    This function calculates the coagulation rate.
          Inputs:
              N = sample size
              particles = matrix to track soot particles
          Outputs:
              The rate of coagulation (1/(m^3*s))
    """
    kernel_1 = 0
    kernel_2a = 0
    kernel_2b = 0
    number_of_particles = 0
    for row in range(0, len(particles)):
        for index in range(0, len(particles[row])):
            if particles[row][index] != 0:
                number_of_particles += particles[row][index]
                kernel_1 += particles[row][index]*(index + 32)**(1./6.)
                kernel_2a += particles[row][index]*(index + 32)**(2./3.)
                kernel_2b += particles[row][index]*(index + 32)**(-1./2.)

    if number_of_particles <= 2:
        rate_of_coagulation = 0
    else:
        rate_of_coagulation = (1.4178/N)*((number_of_particles-2)*kernel_1 + kernel_2a*kernel_2b)
    return rate_of_coagulation


# =========================================== #
# Wait an exponentially distributed time step #
# =========================================== #


def time_step(inception_rate, coagulation_rate):
    """
    This function calculates the exponentially distributed time step.
        Inputs:
            inception_rate = (1/(m^3*s))
            coagulation_rate = (1/(m^3*s))
        Output:
            tao = the exponentially distributed time step
    """
    r = random.uniform(0, 1)
    waiting_parameter = inception_rate + coagulation_rate
    tao = waiting_parameter**(-1)*math.log(1/r)
    return tao


# ============================== #
# Select event probabilistically #
# ============================== #


def select_event(inception_rate, coagulation_rate):
    """
    This function selects the next event probabilistically.
            Inputs:
                inception_rate (1/(m^3*s))
                coagulation_rate (1/(m^3*s))
            Outputs:
                r = 1 indicates inception selection
                r = 2 indicates coagulation selection
    """
    r = random.uniform(0, 1)
    if r >= 0 and r <= (inception_rate/(inception_rate + coagulation_rate)):
        selection = 1
    else:
        selection = 2
    return selection


# ============== #
# Inception Step #
# ============== #


def inception_step(particles):
    """
    This function selects a group probabilistically, then adds a particle of size 32 to the group"
        Inputs:
        Outputs:
    """

    # Determine probabilities for each "group" (e.g. row) of particles
    number_of_particles = 0
    row_sum = 0
    row_sum_list = np.zeros(len(particles))
    for row in range(0, particles-1):
        for index in range(0, particles[i]-1):
            if index != 0:
                number_of_particles += particles[row][index]
                row_sum += particles[row][index]*[index+32]**(1/6)
        row_sum_list[row] = (row_sum)

    total_sum = np.sum(row_sum_list)

    row_probability = np.zeros(len(particles))
    for index in row_sum_list:
        row_probability[index] = row_sum_list[index]/total_sum

    # select group stochastically
    r = random.uniform(0, 1)
    for index in row_probability:
        if r >= sum(row_probability[0:index+1])
            continue
        else:
            group_selection = index - 1
            particles[index - 1][0] += 1
    return print ("Inception Step Complete")


"""def coagulation_step()
"""

# Run Simulation
sample_times = sample_times(start_time, stop_time)
tao = time_step(inception_rate(1900, 0.3), coagulation_rate(N, particles))
event = select_event(inception_rate(1900, 0.3), coagulation_rate(N, particles))

if event = 1:
    inception_step = inception_rate(particles)
