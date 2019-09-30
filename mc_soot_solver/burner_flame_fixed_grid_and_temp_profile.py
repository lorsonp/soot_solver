"""
A burner-stabilized premixed argon-oxygen-acetylene flame at low pressure.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import operator
from get_temp_profile import t_distance_vals, temp_vals, get_temp_prof
from get_conc_profile import *
from temp_and_conc_comparison import distance_grid, fitted_temp, fitted_concentration, c_distance_vals, concentration_vals

# ---------------- #
# Flame Simulation #
# ---------------- #

p = 0.12 * ct.one_atm
tburner = 300.0
reactants = 'AR:0.55, O2:0.2143, C2H2:0.2357, A4:0'  # premixed gas composition
grid = np.linspace(0, 0.09, 20) # m
grid_points = len(grid)
loglevel = 1  # amount of diagnostic output (0 to 5)

# Create gas object
gas = ct.Solution('abf_mech90torr.cti')
gas.TPX = tburner, p, reactants

mdot = 0.204 * gas.density  # kg/(m^2*s)

# Create flame object
f = ct.BurnerFlame(gas, grid=grid)
f.burner.mdot = mdot
f.energy_enabled = False
zloc, tvalues = get_temp_prof(t_distance_vals, temp_vals, points=grid_points)
f.flame.set_fixed_temp_profile(zloc, tvalues)
conc_grid, cvalues = get_conc_prof(c_distance_vals, concentration_vals, points=grid_points)
f.X[gas.species_index('A4'), :] = cvalues
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Mix'
f.solve(loglevel, refine_grid=False)
f.save('c2h2_o2_ar_burner_flame.xml', 'mix', 'solution with mixture-averaged transport')

f.transport_model = 'Multi'
f.solve(loglevel) # don't use 'auto' on subsequent solves
f.show_solution()
f.save('c2h2_o2_ar_burner_flame.xml', 'multi', 'solution with multicomponent transport')

f.write_csv('c2h2_o2_ar_burner_flame.csv', quiet=False)

# Get concentration indices
iA4 = gas.species_index('A4')
iC2H2 = gas.species_index('C2H2')
iCO = gas.species_index('CO')
iH = gas.species_index('H')
iH2 = gas.species_index('H2')
iH2O = gas.species_index('H2O')
iO2 = gas.species_index('O2')
iOH = gas.species_index('OH')

# Did it work? Nope.

# Plot measured/calculated temps and concentrations
fig, (temp_comp, conc_comp) = plt.subplots(1, 2)
fig.suptitle('Measured vs. Calculated Data Comparison (Fixed Temp Profile)')

# Plot temp
temp_comp.plot(zloc, tvalues, 'o', color='g')
temp_comp.plot(distance_grid, fitted_temp, '-', color='g', label='measured curve')
plot_grid = [element/f.grid[-1] for element in f.grid]
temp_comp.plot(plot_grid, f.T, '-', color='b', label='calculated curve')
temp_comp.title.set_text('Temperature (K)')
temp_comp.legend()
plt.show()
