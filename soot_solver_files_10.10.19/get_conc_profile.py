"""This script takes a measured concentration profile of pyrene (published by H. Bockhorn, F. Fetting,and H.-W. Wenz,
1983) and outputs a fixed flame profile, which will then be input into a Cantera 1D burner flame.

Credit for the published temp and pyrene concentration data:
H. Bockhorn, F. Fetting, H.-W. Wenz
Ber. Bunsenges. Phys. Chem., 87 (1983), p. 1067"""

import numpy as np
import math
import matplotlib.pyplot as plt
import csv
from scipy import interpolate


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


def get_conc_prof(c_distance_vals, concentration_vals, width, points):

    # Create distance grid
    distance_grid = np.linspace(0, 1, points)

    # Define ranges for fits
    range_1 = int(round(1000 * (c_distance_vals[12] / width)))
    range_2 = int(round(1000 * ((c_distance_vals[14] - c_distance_vals[12]) / width)))
    range_3 = int(round(1000 * ((c_distance_vals[17] - c_distance_vals[14]) / width)))
    range_4 = int(round(1000 * ((c_distance_vals[-3] - c_distance_vals[17]) / width)))
    range_5 = int(round(1000 * ((c_distance_vals[-1] - c_distance_vals[-3]) / width)))

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
    fitted_concentration = []  # molefraction

    for index in distance_grid:
        if 0 <= index*width < x_vals_1[-1]:
            fitted_concentration.append(10**f_1(index*width))

        elif x_vals_2[0] <= index*width < x_vals_2[-1]:
            fitted_concentration.append(10**interpolate.splev(index*width, splines_2))

        elif x_vals_3[0] <= index*width < x_vals_3[-1]:
            fitted_concentration.append(10**f_3(index*width))

        elif x_vals_4[0] <= index*width < x_vals_4[-1]:
            fitted_concentration.append(10**f_4(index*width))

        elif x_vals_5[0] <= index*width <= 0.09:
            fitted_concentration.append(10**interpolate.splev(index*width, splines_5))

    conc_grid = distance_grid  # [0,1]

    return conc_grid, fitted_concentration

# Get measured Conc data points
d_to_convert = [0.570, 0.577, 0.586, 0.597, 0.608, 0.618, 0.624, 0.633, 0.647, 0.659, 0.673, 0.686, 0.707, 0.735, 0.761,
            0.775, 0.784, 0.809, 0.854, 0.940, 1.080, 1.106, 1.295, 1.484, 1.674, 1.864, 2.054, 2.246]  # inches
c_to_convert = [5.067, 5.0114, 4.9357, 4.860, 4.7837, 4.7074, 4.6311, 4.5548, 4.4785, 4.4022, 4.3259, 4.2496, 4.1733,
                4.136, 4.1733, 4.2496, 4.3259, 4.4022, 4.4785, 4.5548, 4.594, 4.593, 4.584, 4.550, 4.507, 4.454, 4.392, 4.326]  # inches

c_distance_vals, concentration_vals = convert_d_and_c(d_to_convert, c_to_convert)

# conc_grid, cvalues = get_conc_prof(c_distance_vals, concentration_vals, width=0.09, points=1000)

