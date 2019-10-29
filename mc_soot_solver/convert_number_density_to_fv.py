import numpy as np
import csv

MassSoot = 12.0E-3_WP        ! kg/mol
rho_soot = 1800.0_WP         ! kg/m^3
avogadro  = 6.02214199E23_WP ! 1/mol

# --------------------- #
# Calculate Soot Volume #
# --------------------- #

# convert # carbons to diameter
Carbon_2_diam = ((6*MassSoot*number_carbon)/(np.pi*rho_soot*avogadro))**(1/3)   # m
diam_to_radius = Carbon_2_diam/2                                                # m
radius_to_volume = (4/3)*np.pi*diam_to_radius**3                                # m^3
volume = part_num/predict_density
fv  = M100 * MassSoot/(rho_soot*avogadro*volume)
