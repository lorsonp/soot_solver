"""This is a Stochastic Monte Carlo solver for the particle size distribution of soot.
   For this attempt:
    min particle size = 1 (equivalent to 32 Carbon atoms, or 1 monomer)
    max particle size = 10^6 (equivalent to 3.2 x 10^7 Carbon or 10^6 monomers)
    particles are grouped in re-sizeable arrays. Indices are added or deleted depending on the particle event.
    """

import math
import numpy as np
import random
from burner_flame import pyrene_number_concentration, reaction_temp
import matplotlib.pyplot as plt

# to save time while building, temporarily define concentration and temp
pyrene_concentration = pyrene_number_concentration  # number of pyrene molecules/m^3
count_inception_steps = 0
count_coagulation_steps = 0

def initiate_system():
    """This function initiates the particle state and defines the starting parameters.
            Each iteration of the main function waits an exponentially distributed time step, selects an event, and updates the particle system.
            The function cycles through the aforementioned steps until the stop_time is reached
            Inputs:
            Outputs:
        """
    # ================= #
    # PARTICLE TRACKING #
    # ================= #
    # - each array is a separate 'group' with a total of 21 groups
    # - index values denote the particle size (denoted by the # of monomers, 1 monomer = 32 Carbon atoms)
    # - each group holds a specific range of sizes

    # =================== #
    # GROUP ORGANIZATION: #
    # =================== #
    # Group:    Size:
    # 1         0 - 1
    # 2         1 - 2
    # 3         3 - 4
    # 4         5 - 8
    # 5         9 - 16
    # 6         17 - 32
    # 7         33 - 64
    # 8         65 - 128
    # 9         129 - 256
    # 10        257 - 512
    # 11        513 - 1,024
    # 12        1,025 - 2,048
    # 13        2,049 - 4,096
    # 14        4,097 - 8,192
    # 15        8,193 - 16,384
    # 16        16,385 - 32,768
    # 17        32,769 - 65,536
    # 18        65,537 - 131,072
    # 19        131,073 - 262,144
    # 20        262,145 - 524,288
    # 21        524,289 - 1,048,576

    # set initial parameters
    parameters = dict()
    parameters["N"] = 1000  # sample size
    parameters["B"] = 2  # acceptance probability for groups
    parameters["Number of Particles"] = parameters["N"]

    max_size = 10**6  # maximum particle size
    number_of_groups = math.ceil(np.log(max_size)/np.log(parameters["B"]) + 1)

    # define groups of particles
    particles = dict()
    group_size_bins = dict()

    # begin each group as an empty list and define group size bin limits
    for group in range(1, number_of_groups + 1):
        if group == 1:
            group_size_bins[1] = [0, 1]
            particles[group] = []
        else:
            group_size_bins[group] = [parameters["B"] ** (group - 2), parameters["B"] ** (group - 1)]
            particles[group] = []

    # define initial particle distribution
    for i in range(1, parameters["N"] + 1):
        particles[1].append(1)

    # define start, running, and stop time
    start_time = 0.0
    parameters["Running Time"] = start_time
    parameters["Stop Time"] = 1e-19
    time_step = 0

    parameters["Count Time Steps"] = 0
    parameters["Count Inception Steps"] = 0
    parameters["Count Coagulation Steps"] = 0

    return parameters, particles, group_size_bins, time_step


def main(pyrene_concentration, reaction_temp, parameters, particles, group_size_bins, time_step):
    """This is the main function.
        Each iteration of the main function waits an exponentially distributed time step, selects an event, and updates the particle system.
        The function cycles through the aforementioned steps until the stop_time is reached
        Inputs:
        Outputs:
    """
    while parameters["Running Time"] < parameters["Stop Time"]:
        # print("Time:", parameters["Running Time"], "seconds.")
        # ============================================ #
        # Wait an Exponentially Distributed Time Step #
        # ============================================ #
        inception_rate = get_inception_rate(reaction_temp, pyrene_concentration)
        coagulation_rate = get_coagulation_rate(particles, parameters)
        time_step = calculate_time_step(inception_rate, coagulation_rate)
        parameters["Running Time"] += time_step

        # ============================== #
        # Select event probabilistically #
        # ============================== #
        selection = select_event(inception_rate, coagulation_rate)

        # ============================================ #
        # complete event step & update particle system #
        # ============================================ #
        if selection == 1:
            particles, parameters, pyrene_concentration = inception_step(particles, parameters, pyrene_concentration)
        else:
            particles, parameters = coagulation_step(particles, parameters, group_size_bins)

        parameters["Count Time Steps"] += 1
        if parameters["Count Time Steps"] is 1000:
            print("1000 time steps. Running Time: "+parameters["Running Time"]+".")
            parameters["Count Time Steps"] -= 1000

    return particles, parameters


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


def get_coagulation_rate(particles, parameters):
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

    for group in particles:
        if sum(particles[group]) != 0:
            for index in range(0, len(particles[group])):
                kernel_1 += particles[group][index] ** (1. / 6.)
                kernel_2a += particles[group][index] ** (2. / 3.)
                kernel_2b += particles[group][index] ** (-1. / 2.)

    if parameters["Number of Particles"] <= 2:
        coagulation_rate = 0
    else:
        coagulation_rate = (1.4178/parameters["N"])*((parameters["Number of Particles"]-2)*kernel_1 + kernel_2a*kernel_2b)
    return coagulation_rate


def calculate_time_step(inception_rate, coagulation_rate):
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
    tao = waiting_parameter**(-1)*np.log(1/r)
    """
    for times in sample_times:
        if tao <= times:
            continue
        else:
            discrete_waiting_time = times
    """

    return tao


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
    if 0 <= r <= (inception_rate/(inception_rate + coagulation_rate)):
        selection = 1
    else:
        selection = 2
    return selection


def select_kernel(kernel_1, kernel_2):
    r = random.uniform(0, 1)
    if 0 <= r <= (kernel_1 / (kernel_1 + kernel_2)):
        selection = 1
    else:
        selection = 2
    return selection


def select_group_and_index(particles, parameters, kernel_group_sum, E):
    """
    This function calculates the group probabilities and selects a group stochastically
        Inputs:
        Outputs:
    """
    # Determine probabilities for each group
    group_probability = []
    for group in range(len(kernel_group_sum)):
        group_probability.append(kernel_group_sum[group] / sum(kernel_group_sum))

    # select group stochastically
    r = random.uniform(0, 1)
    for group in range(len(group_probability)):
        if r < sum(group_probability[0:group + 1]):
            group_selection = group + 1
            break

    # Select index stochastically
    index_selection = np.random.randint(len(particles[group_selection]))

    # Verify selection acceptance
    r = random.uniform(0, 1)
    while r > particles[group_selection][index_selection] ** E / (parameters["B"] ** (group_selection - 1)) ** E:
        # resample and verify again
        index_selection = np.random.randint(len(particles[group_selection]))
        r = random.uniform(0, 1)

    return group_selection, index_selection


def select_k1_index_2(parameters):
    index_2_selection = np.random.randint(parameters["Number of Particles"])
    group_sum = []
    for group in particles:
        if len(particles[group]) is not 0:
            group_sum.append(len(particles[group]))
            if index_2_selection <= sum(group_sum):
                group_2 = group
                break

    if group_2 > 1:
        index_2 = index_2_selection - sum(group_sum[0:group]) - 1
    else:
        index_2 = index_2_selection - 1
    return group_2, index_2


def inception_step(particles, parameters, pyrene_concentration):
    """
    This function selects a group probabilistically, then adds a particle of size 32 to the group"
        Inputs:
        Outputs:
    """

    # Add particle to group 1
    particles[1].append(1)
    parameters["Number of Particles"] += 1
    pyrene_concentration = pyrene_concentration - 2*parameters["N"]
    parameters["Count Inception Steps"] += 1
    # print("Inception Step Complete")

    return particles, parameters, pyrene_concentration


def coagulation_step(particles, parameters, group_size_bins):
    # Determine probabilities for each kernel and each group
    kernel_1_group_sum = []
    kernel_2_group_sum = []

    for group in particles:
        kernel_1_counter = 0
        kernel_2_counter = 0
        if sum(particles[group]) is not 0:
            for index_1 in range(len(particles[group]) - 1):
                for index_2 in range(index_1 + 1, len(particles[group])):
                    if index_1 is not particles[group][-1]:
                        kernel_1_counter += particles[group][index_1] ** (1 / 6)
                        kernel_2_counter += particles[group][index_1] ** (2 / 3) * particles[group][index_2] ** (-1 / 2)
                    else:
                        kernel_1_counter += particles[group][index_1] ** (1 / 6)

        kernel_1_group_sum.append(kernel_1_counter)
        kernel_2_group_sum.append(kernel_2_counter)

    kernel_1 = 1.4178 / (2 * parameters["N"]) * sum(kernel_1_group_sum)
    kernel_2 = 1.4178 / (2 * parameters["N"]) * sum(kernel_2_group_sum)

    # Select kernel stochastically
    kernel_selection = select_kernel(kernel_1, kernel_2)

    # Select particles
    if kernel_selection == 1:
        # Select group and index 1
        group_1, index_1 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_1_group_sum, E=1 / 6)
        group_2, index_2 = select_k1_index_2(parameters)

        # if group and index selections are identical, resample
        while group_1 == group_2 and index_1 == index_2:
            group_1, index_1 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_1_group_sum, E=1 / 6)
            group_2, index_2 = select_k1_index_2(parameters)

    else:
        # Select group and index 1
        group_1, index_1 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_2_group_sum, E=2 / 3)
        group_2, index_2 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_2_group_sum, E=-1 / 2)

        # if group and index selections are identical, resample
        while group_1 == group_2 and index_1 == index_2:
            group_1, index_1 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_2_group_sum, E=2 / 3)
            group_2, index_2 = select_group_and_index(particles, parameters, kernel_group_sum=kernel_2_group_sum, E=-1 / 2)

    # Determine whether fictitious & update particle system
    size_1 = particles[group_1][index_1]
    size_2 = particles[group_2][index_2]

    coag_kernel = (1 / size_1 + 1 / size_2) ** 0.5 * ((1 / size_1) ** (1 / 3) + (1 / size_2) ** (1 / 3)) ** 2
    maj_kernel = 1.4178 * (size_1 ** -0.5 + size_2 ** -0.5) * (size_1 ** (2 / 3) + size_2 ** (2 / 3))

    r = random.uniform(0, 1)

    if r <= coag_kernel / maj_kernel:

        # Remove 2 selected particles
        del (particles[group_1][index_1])
        del (particles[group_2][index_2])

        # Add new particle to system
        new_particle_size = size_1 + size_2
        for group in range(1, len(group_size_bins) + 1):
            if group_size_bins[group][0] < new_particle_size <= group_size_bins[group][1]:
                particles[group].append(new_particle_size)
                break

        # Update number of particles
        parameters["Number of Particles"] = parameters["Number of Particles"] - 1
        parameters["Count Coagulation Steps"] += 1
        # print("Coagulation Step Complete")

    else:
        print("Coagulation Step Fictitious")

    return particles, parameters


def plot_size_distribution(particles, parameters):
    # Count number of particles of each size
    particle_size = []
    particle_count = []
    for group in particles:
        for index in particles[group]:
            if index not in particle_size:
                particle_size.append(index)
                particle_count.append(particles[group].count(index))

    particle_size_and_count = {key: value for key, value in zip(particle_size, particle_count)}

    # Sort dictionary by particle size
    sorted_particle_size_and_count = {}
    for key in sorted(particle_size_and_count.keys()):
        sorted_particle_size_and_count[key] = particle_size_and_count[key]
    print(sorted_particle_size_and_count)

    # Calculate number density of each particle size
    number_density = {}
    for size in sorted_particle_size_and_count:
        number_density[size] = sorted_particle_size_and_count[size]/parameters["Number of Particles"]

    # Plot size distribution
    plt.plot(list(number_density.keys()), list(number_density.values()), 'bo')
    plt.title("Particle size distribution at %.2E seconds" % parameters["Running Time"])
    plt.xlabel("Particle size (# of monomers)")
    plt.ylabel("Number concentration")
    plt.xscale('log')
    plt.xlim(1, 1e6)
    plt.show()


# Run Simulation
[parameters, particles, group_size_bins, time_step] = initiate_system()
main(pyrene_concentration, reaction_temp, parameters, particles, group_size_bins, time_step)
plot_size_distribution(particles, parameters)

for group in particles:
    print("Group "+str(group)+" has "+str(len(particles[group]))+" particles.")
