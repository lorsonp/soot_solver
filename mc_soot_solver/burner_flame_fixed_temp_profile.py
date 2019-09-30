"""
A burner-stabilized premixed argon-oxygen-acetylene flame at low pressure with a fixed temperature profile.

Credit for published temp and pyrene concentration data:
H. Bockhorn, F. Fetting, H.-W. Wenz
Ber. Bunsenges. Phys. Chem., 87 (1983), p. 1067
"""

import cantera as ct
import math
from get_temp_profile import *
from get_conc_profile import *
import numpy as np
import matplotlib.pyplot as plt
import operator

# ---------------- #
# Flame Simulation #
# ---------------- #
p = 0.12 * ct.one_atm
tburner = 300.0
reactants = 'AR:0.55, O2:0.2143, C2H2:0.2357'  # premixed gas composition
width = 0.09 # m
loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = True

# Create gas object
gas = ct.Solution('abf_mech90torr.cti')
gas.TPX = tburner, p, reactants

mdot = 0.204 * gas.density  # kg/(m^2*s)

# Create flame object
f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.energy_enabled = False
zloc, tvalues = get_temp_prof(t_distance_vals, temp_vals, width=width, points=100)
f.flame.set_fixed_temp_profile(zloc, tvalues)
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Mix'
f.solve(loglevel, refine_grid)
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

# --------------------------------- #
# Temp and Concentration Comparison #
# --------------------------------- #

# Plot measured/calculated temps and concentrations
fig, (temp_comp, conc_comp) = plt.subplots(1, 2)
fig.suptitle('Measured vs. Calculated Temperature and Concentration (Fixed Temp Profile)')

# Plot temp
fitted_temp_grid, fitted_temp = get_temp_prof(t_distance_vals, temp_vals, width=width, points=1000)
temp_comp.plot([element*width for element in zloc], tvalues, 'o', color='g')
temp_comp.plot([element*width for element in fitted_temp_grid], fitted_temp, '-', color='g', label='Measured (Bockhorn et al.)')
temp_comp.plot(f.grid, f.T, '-', color='b', label='Fixed Profile (Cantera burner flame)')
temp_comp.set_xlabel('Distance (m)')
temp_comp.set_ylabel('Temperature (K)')
temp_comp.legend()

# Plot pyrene concentration
fitted_conc_grid, fitted_conc = get_conc_prof(c_distance_vals, concentration_vals, width=width, points=1000)
conc_comp.plot(c_distance_vals, concentration_vals, 'o', color='g')
conc_comp.plot([element*width for element in fitted_conc_grid], [math.log(element, 10) for element in fitted_conc], '-', color='g', label='Measured (Bockhorn et al.)')
conc_comp.plot(f.grid, [math.log(element, 10) for element in f.X[gas.species_index('A4')]], '-', color='b', label='Calculated (Cantera burner flame)')
conc_comp.set_xlabel('Distance (m)')
conc_comp.set_ylabel('Concentration (mole fraction, log scale)')
conc_comp.legend()
plt.show()

# Get pyrene concentration max/min
max_index, max_pyrene_molar_concentration = max(enumerate(f.concentrations[iA4, :]), key=operator.itemgetter(1))  #kmoles/m^3
min_index, min_pyrene_molar_concentration = min(enumerate(f.concentrations[iA4,:]), key=operator.itemgetter(1))  #kmoles/m^3
reaction_temp = f.T[min_index]

# ---------------------------- #
# Plot Major and Minor Species #
# ---------------------------- #

fig, (major_sp, minor_sp) = plt.subplots(1, 2)
fig.suptitle('Major and Minor Species Concentrations (Fixed Temp Profile)')

# Plot major species
major_sp.plot(f.grid, f.X[gas.species_index('C4H2'), :], 'b', label='C4H2')
major_sp.plot(f.grid, f.X[gas.species_index('C6H2'), :], 'g', label='C6H2')
major_sp.plot(f.grid, f.X[gas.species_index('A1'), :], 'r', label='benzene')
major_sp.plot(f.grid, f.X[gas.species_index('C2H4'), :], 'm', label='C2H4')
major_sp.plot(f.grid, f.X[gas.species_index('A2'), :], 'y', label='naphthalene')
major_sp.set_xticks(np.arange(0, .1, .02))
major_sp.set_yscale('log')
major_sp.set_ylim(10e-8, 10e-2)
major_sp.set_xlabel('Distance (m)')
major_sp.set_ylabel('Concentration (molefraction, log scale)')
major_sp.title.set_text('Major Species')
major_sp.legend()

# Plot minor species
minor_sp.plot(f.grid, f.X[gas.species_index('CO'), :], 'b', label='CO')
minor_sp.plot(f.grid, f.X[gas.species_index('H2'), :], 'g', label='H2')
minor_sp.plot(f.grid, f.X[gas.species_index('C2H2'), :], 'r', label='C2H2')
minor_sp.plot(f.grid, f.X[gas.species_index('H2O'), :], 'm', label='H2O')
minor_sp.plot(f.grid, f.X[gas.species_index('CO2'), :], 'y', label='CO2')
minor_sp.plot(f.grid, f.X[gas.species_index('O2'), :], 'c', label='O2')
minor_sp.plot(f.grid, f.X[gas.species_index('CH4'), :], 'k', label='CH4')
minor_sp.set_xticks(np.arange(0, .1, .02))
minor_sp.set_yscale('log')
minor_sp.set_ylim(10e-4, 0)
minor_sp.set_xlabel('Distance (m)')
minor_sp.set_ylabel('Concentration (molefraction, log scale)')
minor_sp.title.set_text('Minor Species')
minor_sp.legend()
plt.show()
