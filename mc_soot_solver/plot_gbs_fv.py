import numpy as np
import csv
import matplotlib.pyplot as plt


def read_data(file, data):
    f = open(file, 'r')
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)
    data = [r for r in reader]

    return data


# Constants
MassSoot = 12.0E-3                  # kg/mol
rho_soot = 1800                     # kg/m^3
avogadro = 6.02214199E23           # 1/mol


# Read in Guillaumes fv data
g_fv = []       # unit-less
g_dist = []     # m
g_time = []     # ms
[g_fv, g_dist, g_time] = read_data(file='gb_data.csv', data=[g_fv, g_dist, g_time])

# Read in Experimental data
exp_fv = []     # unit-less
exp_dist = []   # m
[exp_dist, exp_fv] = read_data(file='ISF_Flame_2B_experimental_data.csv', data=[exp_dist, exp_fv])

# Read in my data
my_particles_to_convert = []                  # number/cm^3
my_part_num_to_convert = []     # number/cm^3
my_time = []                    # s
my_dist = []                    # m
my_num_carbon_to_convert = []   # number/cm^3
[my_particles_to_convert, my_part_num_to_convert, my_time, my_dist, my_num_carbon_to_convert] = read_data(file='ISF_flame_2b_soot_data.csv', data=[my_particles_to_convert, my_part_num_to_convert, my_time, my_dist, my_num_carbon_to_convert])

# Convert my data from number/cm^3 to number/m^3
my_particles = [index*1e6 for index in my_particles_to_convert]         # number/m^3
my_part_num = [index*1e6 for index in my_part_num_to_convert]           # number/m^3
my_num_carbon = [index*1e6 for index in my_num_carbon_to_convert]       # number/m^3
my_predict_dens = max(my_part_num)                                      # number/m^3

# Calculate volume and fv at each data point
my_volume = []
my_fv = []

for index in range(len(my_part_num)):
    my_volume.append(my_part_num[index]/my_predict_dens)
    my_fv.append(my_num_carbon[index]*MassSoot/(rho_soot*avogadro*my_volume[index]))


# plot my num dens data
plt.plot([element*1e3 for element in my_dist], [element/(1e6)**3 for element in my_part_num])
plt.show()
"""    
# OPTION 2: We'll return to this later
# Convert number of carbon to vol
CarbonToDiam = []
for index in my_num_carbon:
    CarbonToDiam.append(((6*MassSoot*index)/(np.pi*rho_soot*avogadro))**(1/3))  # m


# Plot data
plt.plot(g_dist, g_fv, 'g', label='Guillaume\'s data')
plt.plot(exp_dist, exp_fv, linestyle='--', color='k', label='Experimental Data')
plt.plot(my_dist, my_fv, 'b', label='My data')
plt.xlabel('Distance (m)')
plt.ylabel('Volume Fraction (ppm)')
plt.title("Volume Fraction vs. Distance")
plt.legend()
plt.show()

# Truncate data to distance values less than 0.003 m
exp_dist_trunc = [index for index in exp_dist if index<=0.003]
exp_fv_trunc = exp_fv[0:len(exp_dist_trunc)]
g_dist_trunc = [index for index in g_dist if index<=0.003]
g_fv_trunc = g_fv[0:len(g_dist_trunc)]

# Plot data
plt.plot(g_dist_trunc, g_fv_trunc, 'g', label='Guillaume\'s data')
plt.plot(exp_dist_trunc, exp_fv_trunc, linestyle='--', color='k', label='Experimental Data')
plt.plot(my_dist, my_fv, 'b', label='My data')
plt.xlabel('Distance (m)')
plt.ylabel('Volume Fraction')
plt.title("Volume Fraction vs. Distance")
plt.legend()
plt.show()"""