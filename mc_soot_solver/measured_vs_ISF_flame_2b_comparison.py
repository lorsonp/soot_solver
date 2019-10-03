""" This script compares the temp and pyrene concentration profiles of measured values and calculated from ISF Flame 2b
Credit for the input data is attributed to Guillaume Blanquardt, et al."""

import numpy as np
import math
from get_temp_profile import *
from get_conc_profile import *
import numpy as np
import matplotlib.pyplot as plt

# -------- #
# Get Data #
# -------- #

# read in FlameMaster Data
f = open('ISF_flame_2b_data.csv', 'r')
reader = csv.reader(f, quoting=csv.QUOTE_NONNUMERIC)

[ISF_distance, ISF_time, ISF_temp, ISF_conc] = [r for r in reader]  # m, ms, K, molefraction

# get measured temp and conc prof
measured_conc_grid, measured_conc = get_conc_prof(c_distance_vals, concentration_vals, width=ISF_distance[-1], points=100)
measured_temp_grid, measured_temp = get_temp_prof(t_distance_vals, temp_vals, width=ISF_distance[-1], points=100)

# --------------------------- #
# Plot Temp and Concentration #
# --------------------------- #

# plot ISF data vs. measured data
# plot_results(ISF_distance, measured_temp_grid, measured_temp, measured_conc_grid, measured_conc)
def plot_results(ISF_distance, measured_temp_grid, measured_temp, measured_conc_grid, measured_conc):
    fig, (t, c) = plt.subplots(1,2)
    fig.suptitle('Measured vs. Calculated Temperature and Concentration')
    t.plot(ISF_distance, ISF_temp, 'b', label='Calculated (ISF Flame 2B)')
    t.plot([i*ISF_distance[-1] for i in measured_temp_grid], measured_temp, 'g', label='Measured (Bockhorn et al.)')
    t.set_xlabel('Distance (m)')
    t.set_ylabel('Temperature (K)')
    t.legend()

    c.plot(ISF_distance, ISF_conc, 'b', label='Calculated (ISF Flame 2B)')
    c.plot([i*ISF_distance[-1] for i in measured_conc_grid], measured_conc, 'g', label='Measured (Bockhorn et al.)')
    c.set_xlabel('Distance (m)')
    c.set_ylabel('Concentration (molefraction)')
    c.set_yscale('log')
    c.legend(loc='lower right')
    plt.show()

    return t, c