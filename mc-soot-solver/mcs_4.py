import math
import numpy as np
import random
# from burner_flame import max_pyrene_molar_concentration, reaction_temp
import matplotlib.pyplot as plt

max_pyrene_molar_concentration = 2.117959208612895e-08      # kmole/m^3
reaction_temp = 1916.1938344245214                          # K

def initiate_system(max_pyrene_molar_concentration, reaction_temp):
    # constants
    C = {}
    C["NA"] = 6.0243e23     # mol^-1
    C["kB"] = 1.3807e-16    # m^2*kg/(s^2*K)
    C["mc"] = 12            # amu
    C["rho"] = 1.8e-3       # kg/cm^3

    # set initial parameters
    P = {}
    P["pyrene_number_concentration"] = max_pyrene_molar_concentration * 1000 * C["NA"] / 1e6  # molecules/cm^3
    P["temp"] = reaction_temp
    P["N"] = 1000  # sample size
    P["B"] = 2  # acceptance probability for groupsj
    mc_kg = C["mc"]*1.661e-27   # convert to kg
    kB_cm = C["kB"]*1e4         # convert to cm^2*kg/(s^2*K)
    P["A"] = 2.2 * (3*mc_kg/(4*np.pi*C["rho"]))**(1/6) * (6*kB_cm*P["temp"]/C["rho"])  # cm/s

    max_size = 3.2e7  # maximum particle size

    # define particles
    particles = [0]

    # define start, running, and stop time
    start_time = 0.0
    P["Running Time"] = start_time
    P["Stop Time"] = 1

    P["Count Time Steps"] = 0
    P["Count Inception Steps"] = 0
    P["Count Coagulation Steps"] = 0
    P["Count Fictitious Coagulation Steps"] = 0

    return C, P, particles


def main(C, P, particles):

    # while P["Running Time"] < P["Stop Time"]:
    for i in range(10):
        # ============================================ #
        # Wait an Exponentially Distributed Time Step #
        # ============================================ #
        inception_rate = get_inception_rate(C, P)
        coagulation_rate, kernel_1, kernel_2 = get_coagulation_rate(C, P, particles)
        P["Running Time"] += calculate_time_step(rates=[inception_rate, coagulation_rate])

        # ============================== #
        # Select event probabilistically #
        # ============================== #
        selection = select_probability(rates=[inception_rate, coagulation_rate]) + 1

        # ============================================ #
        # complete event step & update particle system #
        # ============================================ #
        selection = 1
        if selection == 1:
            particles, P = inception_step(particles, P)
        else:
            particles, P = coagulation_step(particles, P, Y1, Y2)

    return particles, P


def get_inception_rate(C, P):
    C["mc"] = C["mc"] * 1.661e-27  # convert amu > kg
    C["kB"] = C["kB"] * 1e4  # convert m^2*kg/(s^2*K) to cm^2*kg/(s^2*K)
    d_PAH = 7.9e-8  # cm

    pyrene_coagulation = 2.2 * (np.pi*C["kB"]*P["temp"]/C["mc"])**0.5 * d_PAH**2  # cm^3/s
    inception_rate = 0.5 * pyrene_coagulation * P["pyrene_number_concentration"]**2 * P["N"]

    return inception_rate


def get_coagulation_rate(C, P, particles):
    number_of_particles = sum(particles)

    if number_of_particles <= 2:
        coagulation_rate = 0
    else:
        kernel_1 = np.zeros(len(particles))
        kernel_2 = np.zeros(len(particles))

        for index in range(len(particles)):
            np.put(kernel_1, index, particles[index] * (index + 32)**(1/6))

            kernel_2b = 0
            for index_2 in range(len(particles)):
                kernel_2b += particles[index_2] * (index_2 + 32) ** (-1 / 2)

            np.put(kernel_2, index, particles[index] * (index + 32) ** (2 / 3) \
                   * kernel_2b - particles[index] * (index + 32) ** (1 / 6))

        kernel_sum = (number_of_particles - 1)*sum(kernel_1) + sum(kernel_2)
        coagulation_rate = 1.4178 * P["A"]/P["N"] * kernel_sum

    return coagulation_rate, kernel_1, kernel_2


def calculate_time_step(rates):

    r = random.uniform(0, 1)
    waiting_parameter = sum(rates)

    tao = -np.log(r)/waiting_parameter

    return tao


def select_probability(rates):

    probability = [rate / sum(rates) for rate in rates]

    r = random.uniform(0, 1)
    selection = 0

    for index in range(len(probability)):
        if r <= sum(probability[0:index+1]):
            break
        else:
            selection += 1

    return selection


def inception_step(particles, P):

    # Add particle of size 32 to ensemble
    particles[0] += 1
    P["pyrene_number_concentration"] = P["pyrene_number_concentration"] - 2*P["N"]  # how to determine # of particles per stochastic particle?
    P["Count Inception Steps"] += 1

    return particles, P


def select_particle(particles):

    i_index = random.uniform(1, sum(particles))

    for index in range(len(particles)):
        if i_index <= sum(particles[0:index + 1]):
            i_size = index + 32

    return i_size, i_index


def select_particle_weighted(particles, rates):

    selection = select_probability(rates)
    j_size = selection + 32
    r = np.random.randint(1, particles[selection])
    j_index = sum(particles[0:selection]) + r

    return j_size, j_index


def coagulation_step(particles, P, kernel_1, kernel_2):

    # Select kernel
    kernel_selection = select_probability(rates=[sum(kernel_1), sum(kernel_2)]) + 1

    if kernel_selection == 1:

        # Select particles
        i_size, i_index = select_particle(particles)
        j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E=1/6)

        # Verify i_index != j_index
        while i_index == j_index:
            i_size, i_index = select_particle(particles)
            j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E=1/6)

    else:

        # Select particles
        i_size, i_index = select_particle_weighted(particles, rates=kernel_2, E=1 / 6)
        j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E=-1 / 2)

        # Verify i_index != j_index
        while i_index == j_index:
            i_size, i_index = select_particle_weighted(particles, rates=kernel_2, E=1 / 6)
            j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E=-1 / 2)

    # Determine whether fictitious & update particle system

    coag_kernel = (1 / i_size + 1 / j_size) ** 0.5 * (i_size ** (1 / 3) + j_size ** (1 / 3)) ** 2
    maj_kernel = 1.4178 * (i_size ** -0.5 + j_size ** -0.5) * (i_size ** (2 / 3) + j_size ** (2 / 3))

    r = random.uniform(0, 1)

    if r <= coag_kernel / maj_kernel:

        # Remove 2 selected particles
        particles[i_size - 32] = particles[i_size - 32] - 1
        particles[j_size - 32] = particles[j_size - 32] - 1

        # Add new particle to system
        new_particle_size = i_size + j_size
        if new_particle_size - 32 > len(particles) - 1:
            for index in range(len(particles), new_particle_size - 32):
                particles[index] = 0
            particles[new_particle_size - 32] = 1
        else:
            particles[new_particle_size - 32] += 1

        parameters["Count Coagulation Steps"] += 1

    else:
        parameters["Count Fictitious Coagulation Steps"] += 1

    return particles


[C, P, particles] = initiate_system(max_pyrene_molar_concentration, reaction_temp)
main(C, P, particles)