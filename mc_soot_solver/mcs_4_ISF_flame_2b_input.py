# This version of the solver takes input from ISF Flame 2B.  Credit for the input data is attributed to Guillaume Blanquardt, et al.

import math
import numpy as np
import random
from measured_vs_ISF_flame_2b_comparison import *
import matplotlib.pyplot as plt
import operator
from scipy.interpolate import interp1d

ISF_time_grid = [element/1000 for element in ISF_time]  # convert ms to s
def initiate_system(pyrene_molar_concentration, reaction_temp, time_grid, distance_grid):
    # constants
    C = {}
    C["NA"] = 6.0221409e23          # mol^-1
    C["kB"] = 1.38064852e-23 * 1e4  # convert m^2*kg/(s^2*K) to cm^2*kg/(s^2*K)
    C["mc"] = 12 * 1.661e-27        # convert amu > kg
    C["rho"] = 1.8e-3               # kg/cm^3

    # set initial parameters
    P = {}
    P["pyrene_number_concentration"] = [element * 1000 * C["NA"] / 1e6 for element in pyrene_molar_concentration]  # molecules/cm^3
    P["temp"] = reaction_temp
    P["time_grid"] = time_grid
    P["N"] = 1e3  # sample size
    max_size = 3.2e7  # maximum particle size

    # define particles
    particles = [0]

    # define functions to interpolate temp and number concentration
    F = {}
    F["interp_temp_function"] = interp1d(P["time_grid"], P["temp"])
    F["interp_conc_function"] = interp1d(P["time_grid"], P["pyrene_number_concentration"])

    # define start, running, and stop time (s)
    start_time = 0.0
    P["Running Time"] = start_time
    P["Stop Time"] = 0.01  # cannot exceed 0.112 s

    P["Count Time Steps"] = 0
    P["Count Inception Steps"] = 0
    P["Count Coagulation Steps"] = 0
    P["Count Fictitious Coagulation Steps"] = 0

    # define arrays to hold plotting values
    soot_number_density = []  # total number of stochastic particles per unit volume
    time_plot_points = []      # s
    distance_plot_points = []  # m
    tao = []
    count_time = 0
    return C, F, P, particles, tao, count_time, soot_number_density, time_plot_points, distance_plot_points


def main(C, F, P, particles, tao, count_time, soot_number_density, time_plot_points, distance_plot_points):

    while P["Running Time"] < P["Stop Time"]:
    # for i in range(1000000):
        # ================================= #
        # Get temp and conc at current time #
        # ================================= #
        interpolated_temp, interpolated_pyrene_conc = interpolate_value(P, F)

        # ============================================ #
        # Wait an Exponentially Distributed Time Step #
        # ============================================ #
        inception_rate = get_inception_rate(C, P, interpolated_temp, interpolated_pyrene_conc)
        coagulation_rate, kernel_1, kernel_2 = get_coagulation_rate(C, P, interpolated_temp, particles)
        tao = calculate_time_step(rates=[inception_rate, coagulation_rate])
        P["Running Time"] += tao
        count_time += tao

        # ============================== #
        # Select event probabilistically #
        # ============================== #
        selection = select_probability(rates=[inception_rate, coagulation_rate]) + 1

        # ============================================ #
        # complete event step & update particle system #
        # ============================================ #

        if selection == 1:
            particles, P = inception_step(particles, P)
        else:
            particles, P = coagulation_step(particles, P, kernel_1, kernel_2)

        P["Count Time Steps"] += 1
        if P["Count Time Steps"] == 100000:
            # print("time step: "+str(tao))
            # print("pyrene number concentration: "+str(P["pyrene_number_concentration"]))
            print("100,000 time steps. Running time: "+str(P["Running Time"])+" seconds.")
            P["Count Time Steps"] = P["Count Time Steps"] - 10000

        if count_time >= 5e-5:
            velocity = 0.0673  # m/s
            soot_number_density.append(sum(particles))  # total number of stochastic particles per unit volume
            time_plot_points.append(P["Running Time"])  # s
            distance_plot_points.append(P["Running Time"] * velocity)  # m
            count_time = count_time - 5e-5

    print(P["Running Time"])

    return particles, P, soot_number_density, time_plot_points, distance_plot_points

def interpolate_value(P, F):

    interpolated_temp = float(F["interp_temp_function"](P["Running Time"]))
    interpolated_conc = float(F["interp_conc_function"](P["Running Time"]))

    return interpolated_temp, interpolated_conc

def get_inception_rate(C, P, interpolated_temp, interpolated_pyrene_concentration):

    d_PAH = 7.9e-8  # cm
    pyrene_coagulation = 2.2 * (np.pi*C["kB"]*interpolated_temp/C["mc"])**0.5 * d_PAH**2  # cm^3/s
    inception_rate = 0.5 * pyrene_coagulation * interpolated_pyrene_concentration**2

    # Adjust inception_rate so 20 pyrenes make 10 soots
    inception_rate = inception_rate/10

    return inception_rate


def get_coagulation_rate(C, P, interpolated_temp, particles):
    A = 2.2 * (3 * C["mc"] / (4 * np.pi * C["rho"])) ** (1 / 6) * (6 * C["kB"] * interpolated_temp / C["rho"])  # cm/s
    number_of_particles = sum(particles)

    if number_of_particles <= 2:
        coagulation_rate = 0
        kernel_1 = 0
        kernel_2 = 0
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
        coagulation_rate = 1.4178 * A/P["N"] * kernel_sum

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
    particles[0] += 10

    P["Count Inception Steps"] += 1

    return particles, P


def select_particle(particles):

    i_index = random.uniform(1, sum(particles))

    for index in range(len(particles)):
        if i_index <= sum(particles[0:index + 1]):
            i_size = index + 32

    return i_size, i_index


def select_particle_weighted(particles, rates, E):

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
        j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E= -1 / 2)

        # Verify i_index != j_index
        while i_index == j_index:
            i_size, i_index = select_particle_weighted(particles, rates=kernel_2, E=1 / 6)
            j_size, j_index = select_particle_weighted(particles, rates=kernel_1, E= -1 / 2)

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

        P["Count Coagulation Steps"] += 1

    else:
        P["Count Fictitious Coagulation Steps"] += 1

    return particles, P

# def plot_results():

# Plot number density (cm^-3) over time (ms)


# -------------- #
# Run Simulation #
# -------------- #

[C, F, P, particles, tao, count_time, soot_number_density, time_plot_points, distance_plot_points] = initiate_system(pyrene_molar_concentration=ISF_conc, reaction_temp=ISF_temp, time_grid=ISF_time_grid, distance_grid=ISF_distance)
[particles, P, soot_number_density, time_plot_points, distance_plot_points] = main(C, F, P, particles, tao, count_time, soot_number_density, time_plot_points, distance_plot_points)

# ------------ #
# Plot Results #
# ------------ #
plt.plot(time_plot_points, soot_number_density, 'k', label='Number Density')
plt.yscale('log')
plt.show()

print(particles)
print(P)
with open('ISF_flame_2b_soot_data.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerows([particles, soot_number_density, time_plot_points, distance_plot_points])
