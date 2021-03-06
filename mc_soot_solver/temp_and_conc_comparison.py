"""This script takes a measured temperature and concentration profile of pyrene (published by H. Bockhorn, F. Fetting,
and H.-W. Wenz, 1983), fits curves to the data, and compares the fitted measured data to output from a Cantera 1D burner
flame. The array of temps and concentrations with respect to distance can then be used as input for the soot solver.

Credit for the published temp and pyrene concentration data:
H. Bockhorn, F. Fetting, H.-W. Wenz
Ber. Bunsenges. Phys. Chem., 87 (1983), p. 1067"""

import numpy as np
import math
import matplotlib.pyplot as plt
import csv
from scipy import interpolate
# from get_temp_profile import *
# from get_conc_profile import *


def get_data(file):

    # import distance grid, temp, and concentration
    f = open(file, 'r')
    reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)

    [distance_grid, calculated_temp, calculated_concentration] = [r for r in reader]

    new_concentration = []

    for element in calculated_concentration:
        new_concentration.append(math.log(element, 10))

        # calculated_concentration = [math.log(element, 10) for element in calculated_concentration]

    return distance_grid, calculated_temp, calculated_concentration, new_concentration


def get_fitted_temp(t_distance_vals, temp_vals, distance_grid):

    # Define ranges for fits
    range_1 = round(1000*(t_distance_vals[1]/.09))
    range_2 = round(1000*(t_distance_vals[2]-t_distance_vals[1])/.09)
    range_3 = round(1000*((.09-t_distance_vals[2])/.09))

    # Create fits
    x_vals = np.linspace(0, .09, 1000)
    x_vals_1 = np.linspace(0, t_distance_vals[1], range_1)
    x_vals_2 = np.linspace(t_distance_vals[1], t_distance_vals[2], range_2)
    x_vals_3 = np.linspace(t_distance_vals[2], .09, range_3)

    x_1 = np.array([t_distance_vals[0], t_distance_vals[1]])
    y_1 = np.array([temp_vals[0], temp_vals[1]])

    x_2 = np.array([t_distance_vals[1], t_distance_vals[2]])
    y_2 = np.array([temp_vals[1], temp_vals[2]])

    f_1 = np.poly1d(np.polyfit(x_1, y_1, 1))
    f_2 = np.poly1d(np.polyfit(x_2, y_2, 2))
    f_3 = np.poly1d(np.polyfit(t_distance_vals[2:-1], temp_vals[2:-1], 2))

    # interpolate values
    fitted_temp = []
    for index in distance_grid:
        if 0 <= index < t_distance_vals[1]:
            fitted_temp.append(f_1(index))
        elif t_distance_vals[1] <= index < t_distance_vals[2]:
            fitted_temp.append(f_2(index))
        elif t_distance_vals[2] <= index:
            fitted_temp.append(f_3(index))

    return fitted_temp


def convert_d_and_c(d_to_convert, c_to_convert):
    c_distance_vals = np.zeros(len(d_to_convert))
    concentration_vals = np.zeros(len(c_to_convert))

    # convert distance to mm
    for index in range(len(d_to_convert)):
        c_distance_vals[index] = ((d_to_convert[index] - 0.539) / 18.967)  # m

    c_distance_vals[-1] = 0.09

    # convert concentration to molefraction
    for index in range(len(c_to_convert)):
        if c_to_convert[index] >= 4.860:
            concentration_vals[index] = 10 ** -((c_to_convert[index] - 4.860) / 0.757) - 7  # log power of molefraction
        else:
            concentration_vals[index] = 10 ** -((c_to_convert[index] - 4.097) / 0.763) - 6  # log power of molefraction

    return c_distance_vals, concentration_vals


def get_fitted_concentration(c_distance_vals, concentration_vals, distance_grid):

    # Define ranges for fits
    range_1 = int(round(1000 * (c_distance_vals[12] / .09)))
    range_2 = int(round(1000 * ((c_distance_vals[14] - c_distance_vals[12]) / .09)))
    range_3 = int(round(1000 * ((c_distance_vals[17] - c_distance_vals[14]) / .09)))
    range_4 = int(round(1000 * ((c_distance_vals[-3] - c_distance_vals[17]) / .09)))
    range_5 = int(round(1000 * ((c_distance_vals[-1] - c_distance_vals[-3]) / .09)))

    x_vals_1 = np.linspace(0, c_distance_vals[12], range_1)
    x_vals_2 = np.linspace(c_distance_vals[12], c_distance_vals[14], range_2)
    x_vals_3 = np.linspace(c_distance_vals[14], c_distance_vals[17], range_3)
    x_vals_4 = np.linspace(c_distance_vals[17], c_distance_vals[-3], range_4)
    x_vals_5 = np.linspace(c_distance_vals[-3], c_distance_vals[-1], range_5)

    # Create fits
    f_1 = np.poly1d(np.polyfit(c_distance_vals[0:13], concentration_vals[0:13], 3))
    splines_2 = interpolate.splrep(c_distance_vals, concentration_vals)
    f_3 = np.poly1d(np.polyfit(c_distance_vals[14:19], concentration_vals[14:19], 3))
    f_4 = np.poly1d(np.polyfit(c_distance_vals[17:-1], concentration_vals[17:-1], 6))
    splines_5 = interpolate.splrep(c_distance_vals, concentration_vals)

    # interpolate values
    fitted_concentration = np.zeros(len(distance_grid))

    for index in range(len(distance_grid)):
        if x_vals_1[0] <= distance_grid[index] < x_vals_1[-1]:
            fitted_concentration[index] = f_1(distance_grid[index])

        elif x_vals_2[0] <= distance_grid[index] < x_vals_2[-1]:
            fitted_concentration[index] = interpolate.splev(distance_grid[index], splines_2)

        elif x_vals_3[0] <= distance_grid[index] < x_vals_3[-1]:
            fitted_concentration[index] = (f_3(distance_grid[index]))

        elif x_vals_4[0] <= distance_grid[index] < x_vals_4[-1]:
            fitted_concentration[index] = (f_4(distance_grid[index]))

        elif x_vals_5[0] <= distance_grid[index] <= 0.09:
            fitted_concentration[index] = interpolate.splev(distance_grid[index], splines_5)

    return fitted_concentration


# Get data from cantera
distance_grid, calculated_temp, calculated_concentration, new_concentration = get_data(file='burner_flame_output_data.csv')  # m, K, molefraction

# Measured Temp data points
t_distance_vals = [0.0, 0.00381, 0.0063, 0.01286, 0.02053, 0.02921, 0.03852, 0.04815, 0.05889, 0.07074, 0.08296, 0.09]  # m
temp_vals = [300, 1945.52, 1992, 1950, 1900, 1850, 1800, 1750, 1700, 1650, 1600, 1571]  # K

# Get measured Conc data points
d_to_convert = [0.570, 0.577, 0.586, 0.597, 0.608, 0.618, 0.624, 0.633, 0.647, 0.659, 0.673, 0.686, 0.707, 0.735, 0.761,
            0.775, 0.784, 0.809, 0.854, 0.940, 1.080, 1.106, 1.295, 1.484, 1.674, 1.864, 2.054, 2.246]  # inches
c_to_convert = [5.067, 5.0114, 4.9357, 4.860, 4.7837, 4.7074, 4.6311, 4.5548, 4.4785, 4.4022, 4.3259, 4.2496, 4.1733,
                4.136, 4.1733, 4.2496, 4.3259, 4.4022, 4.4785, 4.5548, 4.594, 4.593, 4.584, 4.550, 4.507, 4.454, 4.392, 4.326]  # inches

c_distance_vals, concentration_vals = convert_d_and_c(d_to_convert, c_to_convert)

# Fit Temp data to Cantera Grid
fitted_temp = get_fitted_temp(t_distance_vals, temp_vals, distance_grid)

# Fit Concentration data to Cantera Grid
fitted_concentration = get_fitted_concentration(c_distance_vals, concentration_vals, distance_grid)

# Plot Calculated against measured Data
fig, (temp_comp, conc_comp) = plt.subplots(1, 2)
fig.suptitle('Measured vs. Calculated Temperature and Concentration')

# Plot temp
temp_comp.plot(t_distance_vals, temp_vals, 'o', color='g')
temp_comp.plot(distance_grid, fitted_temp, '-', color='g', label='Measured (Bockhorn et al.)')
temp_comp.plot(distance_grid, calculated_temp, '-', color='b', label='Calculated (Cantera burner flame)')
temp_comp.set_xlabel('Distance (m)')
temp_comp.set_ylabel('Temperature (K)')
temp_comp.legend()

# Plot concentration
conc_comp.plot(c_distance_vals, concentration_vals, 'o', color='g')
conc_comp.plot(distance_grid, fitted_concentration, '-', color='g', label='Measured (Bockhorn et al.)')
conc_comp.plot(distance_grid, new_concentration, '-', color='b', label='Calculated (Cantera burner flame)')
conc_comp.set_xlabel('Distance (m)')
conc_comp.set_ylabel('Concentration (mole fraction, log scale)')
conc_comp.legend()
plt.show()

