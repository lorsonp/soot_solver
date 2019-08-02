"""This is a Stochastic Monte Carlo solver for the particle size distribution of soot."""

import math
import numpy as np
import random
import cantera as ct
import matplotlib.pyplot as plt

def main():
    """This is the main function, which is called to run the simulation
        Inputs:
        Outputs:
    """
    # ======================== #
    # Generate particle system #
    # ======================== #
    # Import source terms:
    # - Need Temperature and pyrene concentration.
    # - If applicable, import initial particle distribution

    gas1 = ct.Solution('gri30.xml')

    # initiate matrix to track particles
    # - each row is a separate 'group'
    # - column index + 32 indicates the size of the particle (# of carbon atoms)
    # - index values denote the number of particles of that size
    # particles = np.zeros((100, 100))
    particles = [
        [4, 1, 0],
        [4, 1, 0]
    ]

    N = 100  # sample size

    # Generate sample times
    start_time = 0.0
    stop_time = 10.0
    sample_times = get_sample_times(start_time, stop_time)
    time_steps = []
    running_time = start_time
    # ============================================ #
    # Wait an Exponentially Distributed Time Step #
    # ============================================ #
    inception_rate = get_inception_rate(1900, 0.3)
    coagulation_rate = get_coagulation_rate(N, particles)
    time_steps.append(calculate_time_step(inception_rate, coagulation_rate, sample_times))
    running_time += sum(time_steps)

    # ============================== #
    # Select event probabilistically #
    # ============================== #
    selection = select_event(inception_rate, coagulation_rate, N, particles)
    selection = 2
    if selection == 1:
        step = inception_step(particles)
    else:
        step = coagulation_step(particles, N)


def get_sample_times(start_time, stop_time):
    sample_times = np.arange(start_time, stop_time, 0.05)
    return sample_times


def get_inception_rate(reaction_temp, pyrene_concentration):
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
    inception_rate = 0.5*pyrene_coagulation*pyrene_concentration**2.
    return inception_rate


def get_coagulation_rate(N, particles):
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
        coagulation_rate = 0
    else:
        coagulation_rate = (1.4178/N)*((number_of_particles-2)*kernel_1 + kernel_2a*kernel_2b)
    return coagulation_rate


def calculate_time_step(inception_rate, coagulation_rate, sample_times):
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
    for times in sample_times:
        if tao <= times:
            continue
        else:
            discrete_waiting_time = times

    return discrete_waiting_time


def select_event(inception_rate, coagulation_rate, N, particles):
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
    if 0 <= r <= (inception_rate/(inception_rate + coagulation_rate)):
        selection = 1
    else:
        selection = 2
    return selection


def group_selection(particles, E):
    """
    This function calculates the group probabilities and selects a group stochastically
        Inputs:
        Outputs:
    """
    # Determine probabilities for each "group" (e.g. row) of particles
    number_of_particles = 0
    group_prob_sum_list = np.zeros(len(particles))
    for group in range(len(particles)):
        group_prob_sum = 0
        for index in range(len(particles[0])):
            if particles[group][index] != 0:
                number_of_particles += particles[group][index]
                group_prob_sum += particles[group][index] * (index + 32)**E
        group_prob_sum_list[group] = group_prob_sum
    total_sum = np.sum(group_prob_sum_list)

    group_probability = np.zeros(len(particles))
    for group in range(len(group_prob_sum_list)):
        group_probability[group] = group_prob_sum_list[group] / total_sum

    # select group stochastically
    r = random.uniform(0, 1)
    for group in range(len(group_probability)):
        if r < sum(group_probability[0:group + 1]):
            group_selection = group
            break
        else:
            continue
    return group_selection


def inception_step(particles):
    """
    This function selects a group probabilistically, then adds a particle of size 32 to the group"
        Inputs:
        Outputs:
    """

    group = group_selection(particles, E=1)
    particles[group][0] += 1
    return print ("Inception Step Complete")


def coagulation_step(particles, N):
    # Determine probabilities for each kernel
    kernel_1_sum = 0
    kernel_2_sum = 0
    number_of_particles = 0
    for group in range(len(particles)):
        for i in range(len(particles[0]) - 1):
            if particles[group][i] != 0:
                number_of_particles += particles[group][i]
                kernel_1_sum += particles[group][i] * (i + 32) ** (1 / 6)
                for j in range(i + 1, len(particles[0])):
                    if i == len(particles[0]) - 2 and j == len(particles[0]) - 1:
                        number_of_particles += particles[group][j]
                    if particles[group][j] != 0:
                        kernel_2_sum += particles[group][i] * (i + 32) ** (2 / 3) * particles[group][i + 1] * (
                                    j + 32) ** (-1 / 2)
    kernel_1 = 1.4178 / (2 * N) * kernel_1_sum
    kernel_2 = 1.4178 / (2 * N) * kernel_2_sum

    # Select kernel stochastically
    r = random.uniform(0, 1)
    if 0 <= r <= (kernel_1 / (kernel_1 + kernel_2)):
        selection = 1
    else:
        selection = 2

    selection = 1
    # Select group stochastically
    if selection == 1:
        group = group_selection(particles, E=1 / 6)

        # Select index stochastically
        B = 2  # upper bound of probability density function
        index = np.random.randint(len(particles[group]))

        # Verify selection acceptance
        r = random.uniform(0, 1)
        if 0 <= r <= particles[group][index] ** (1 / 6) / (B ** (len(particles[group]) - 2)) ** (1 / 6):
            verify = 1
        else:
            verify = 2

        while verify == 2:
            index = np.random.randint(len(particles[group]))

            # Verify selection acceptance
            r = random.uniform(0, 1)
            if 0 <= r <= particles[group][index] ** (1 / 6) / (B ** (len(particles[group]) - 2)) ** (1 / 6):
                verify = 1
            else:
                verify = 2
    return group, index

# Run Simulation
main()
