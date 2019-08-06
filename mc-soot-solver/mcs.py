"""This is a Stochastic Monte Carlo solver for the particle size distribution of soot."""

import math
import numpy as np
import random
import cantera as ct
import matplotlib.pyplot as plt

def initiate_system():
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
    particles = np.zeros([2, 100], dtype = int)

    N = 1000  # sample size

    start_time = 0.0
    stop_time = 5
    sample_times = get_sample_times(start_time, stop_time)
    running_time = start_time
    time_steps = []

    return particles, N, sample_times, running_time, stop_time, time_steps

def main(particles, N, sample_times, running_time, stop_time, time_steps):
    """This is the main function.
        Each iteration of the main function waits an exponentially distributed time step, selects an event, and updates the particle system.
        The function cycles through the aforementioned steps until the stop_time is reached
        Inputs:
        Outputs:
    """
    while running_time < stop_time:
        print("Time:", running_time, "seconds.")
        # ============================================ #
        # Wait an Exponentially Distributed Time Step #
        # ============================================ #
        inception_rate = get_inception_rate(1900, 1e10)
        coagulation_rate = get_coagulation_rate(N, particles)
        time_steps.append(calculate_time_step(inception_rate, coagulation_rate, sample_times))
        running_time += sum(time_steps)

        # ============================== #
        # Select event probabilistically #
        # ============================== #
        selection = select_event(inception_rate, coagulation_rate, N, particles)

        # ============================================ #
        # complete event step & update particle system #
        # ============================================ #
        if selection == 1:
            particles = inception_step(particles)
        else:
            particles = coagulation_step(particles, N)

    return(particles)


def get_sample_times(start_time, stop_time):
    sample_times = np.arange(start_time, stop_time, 1e-5)
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
    pyrene_coagulation = (2.2*4*2/3)*dA**2*(16)*(np.pi*Boltzmann_k*reaction_temp/(2*reduced_m_pyrene))**0.5
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


def select_group(particles, E):
    """
    This function calculates the group probabilities and selects a group stochastically
        Inputs:
        Outputs:
    """
    # Determine probabilities for each "group" (e.g. row) of particles
    group_prob_sum = 0
    group_prob_sum_list = []
    for group in range(len(particles)):
        for index in range(len(particles[group])):
            group_prob_sum += particles[group][index] * (index + 32)**E
        group_prob_sum_list.append(group_prob_sum)
    total_sum = np.sum(group_prob_sum_list)

    group_probability = []
    for group in range(len(group_prob_sum_list)):
        group_probability.append(group_prob_sum_list[group] / total_sum)

    # select group stochastically
    r = random.uniform(0, 1)
    for group in range(len(group_probability)):
        if r < sum(group_probability[0:group + 1]):
            group_selection = group
            break
    return group_selection


def select_kernel(kernel_1, kernel_2):
    r = random.uniform(0, 1)
    if 0 <= r <= (kernel_1 / (kernel_1 + kernel_2)):
        selection = 1
    else:
        selection = 2
    return selection


def select_index(particles, group, E):
    # Select index stochastically
    B = 2  # upper bound of probability density function
    index = np.random.randint(len(particles[group]))

    # Verify selection acceptance
    r = random.uniform(0, 1)
    while r > (index + 1)**E / (B ** (len(particles[group]) - 2)) ** E:
        index = np.random.randint(len(particles[group]))

        # Verify selection again
        r = random.uniform(0, 1)

    return index


def get_kernel_1_index_2(total_number_of_particles, particles, particles_in_group):
    selected_particle = np.random.randint(1, total_number_of_particles + 1)

    for group in range(len(particles)):
        if selected_particle <= sum(particles_in_group[0:group + 1]):
            group_2 = group
            break

    for index in range(len(particles[group])):
        if index <= sum(particles[group][0:index + 1]):
            index_2 = index
    return [group_2, index_2]


def inception_step(particles):
    """
    This function selects a group probabilistically, then adds a particle of size 32 to the group"
        Inputs:
        Outputs:
    """
    # Select group stochastically
    group = np.random.randint(len(particles))

    # Add particle to selected group
    particles[group][0] += 1
    print("Inception Step Complete")
    return particles


def coagulation_step(particles, N):
    # Determine probabilities for each kernel and tally particle numbers
    kernel_1_sum = 0
    kernel_2_sum = 0
    particles_in_group = []
    for group in range(len(particles)):
        particles_in_group.append(sum(particles[group]))
        for i in range(len(particles[group]) - 1):
            kernel_1_sum += particles[group][i] * (i + 32) ** (1 / 6)
            for j in range(i + 1, len(particles[group])):
                if particles[group][j] != 0:
                    kernel_2_sum += particles[group][i] * (i + 32) ** (2 / 3) * particles[group][i + 1] * (
                                j + 32) ** (-1 / 2)
    kernel_1 = 1.4178 / (2 * N) * kernel_1_sum
    kernel_2 = 1.4178 / (2 * N) * kernel_2_sum
    total_number_of_particles = sum(particles_in_group)

    # Select kernel stochastically
    kernel_selection = select_kernel(kernel_1, kernel_2)

    # Select group & indices stochastically
    if kernel_selection == 1:
        # Use kernel 1 to select index_1 stochastically
        group_1 = select_group(particles, E=1/6)
        index_1 = select_index(particles, group=group_1, E=1/6)

        # Select index_2 stochastically from total number of particles
        [group_2, index_2] = get_kernel_1_index_2(total_number_of_particles, particles, particles_in_group)

        while index_1 == index_2:
            # Reselect index_1
            group_1 = select_group(particles, E=1/6)
            index_1 = select_index(particles, group=group_1, E=1/6)

            # Reselect index_2
            [group_2, index_2] = get_kernel_1_index_2(total_number_of_particles, particles, particles_in_group)
    else:
        # Use kernel 2 to select index_1 stochastically
        group_1 = select_group(particles, E=2/3)
        index_1 = select_index(particles, group=group_1, E=2/3)

        # Use kernel 2 to select index_2 stochastically
        group_2 = select_group(particles, E=-1/2)
        index_2 = select_index(particles, group=group_2, E=-1/2)

        while index_1 == index_2:
            # Reselect index_1
            group_1 = select_group(particles, E=2/3)
            index_1 = select_index(particles, group=group_1, E=2/3)

            # Reselect index_2
            group_2 = select_group(particles, E=-1/2)
            index_2 = select_index(particles, group=group_2, E=-1/2)

    # Determine whether fictitious & update particle system
    size_1 = index_1 + 1-
    size_2 = index_2 + 1
    coag_kernel = (1/size_1 + 1/size_2)**0.5 * ((1/size_1)**(1/3) + (1/size_2)**(1/3))**2
    maj_kernel = 1.4178 * (size_1**-0.5 + size_2**-0.5) * (size_1**(2/3) + size_2**(2/3))
    r = random.uniform(0, 1)
    if r <= coag_kernel/maj_kernel:
        # Remove 1 particle from each index
        particles[group_1][index_1] = particles[group_1][index_1] - 1
        particles[group_2][index_2] = particles[group_2][index_2] - 1

        # Add coagulated particle to system
        group = np.random.randint(len(particles))
        particles[group][index_1 + index_2 + 1] += 1
        print("Coagulation Step Complete")
    else:
        print("Coagulation Step Fictitious")
    return particles


# Run Simulation
[particles, N, sample_times, running_time, stop_time, time_steps] = initiate_system()
main(particles, N, sample_times, running_time, stop_time, time_steps)
