"""This script takes a measured temperature profile (published by H. Bockhorn, F. Fetting, and H.-W. Wenz, 1983)
and outputs a fixed flame profile, which will then be input into a Cantera 1D burner flame.

Credit for the published temp and pyrene concentration data:
H. Bockhorn, F. Fetting, H.-W. Wenz
Ber. Bunsenges. Phys. Chem., 87 (1983), p. 1067"""

import numpy as np
import math
import matplotlib.pyplot as plt
import csv
from scipy import interpolate

def get_temp_prof(t_distance_vals, temp_vals, width, points):

    # Create distance grid
    distance_grid = np.linspace(0, 1, points)

    # Define ranges for fits
    range_1 = round(points*(t_distance_vals[1]/width))
    range_2 = round(points*(t_distance_vals[2]-t_distance_vals[1])/width)
    range_3 = round(points*((width-t_distance_vals[2])/width))

    # Create fits
    x_vals = np.linspace(0, width, points)
    x_vals_1 = np.linspace(0, t_distance_vals[1], range_1)
    x_vals_2 = np.linspace(t_distance_vals[1], t_distance_vals[2], range_2)
    x_vals_3 = np.linspace(t_distance_vals[2], width, range_3)

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
        if 0 <= index*width < t_distance_vals[1]:
            fitted_temp.append(f_1(index*width))
        elif t_distance_vals[1] <= index*width < t_distance_vals[2]:
            fitted_temp.append(f_2(index*width))
        elif t_distance_vals[2] <= index*width:
            fitted_temp.append(f_3(index*width))

    return distance_grid, fitted_temp

# Measured Temp data points
t_distance_vals = [0.0, 0.00381, 0.0063, 0.01286, 0.02053, 0.02921, 0.03852, 0.04815, 0.05889, 0.07074, 0.08296, 0.09]  # m
temp_vals = [300, 1945.52, 1992, 1950, 1900, 1850, 1800, 1750, 1700, 1650, 1600, 1571]  # K

# zloc, tvalues = get_temp_prof(t_distance_vals, temp_vals, points=100)
