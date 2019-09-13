"""
A burner-stabilized premixed argon-oxygen-acetylene flame at low pressure with a fixed temperature profile.

Credit for published temp and pyrene concentration data:
H. Bockhorn, F. Fetting, H.-W. Wenz
Ber. Bunsenges. Phys. Chem., 87 (1983), p. 1067
"""

import cantera as ct
import math
from get_temp_profile import zloc, tvalues
from temp_and_conc_comparison import distance_grid, fitted_temp, fitted_concentration, c_distance_vals, concentration_vals
import numpy as np
import matplotlib.pyplot as plt
import operator


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
f.flame.set_fixed_temp_profile(zloc, tvalues)
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Mix'
f.solve(loglevel, refine_grid)
f.save('c2h2_o2_ar_burner_flame.xml', 'mix', 'solution with mixture-averaged transport')

"""f.transport_model = 'Multi'
f.solve(loglevel) # don't use 'auto' on subsequent solves
f.show_solution()
f.save('c2h2_o2_ar_burner_flame.xml', 'multi', 'solution with multicomponent transport')

f.write_csv('c2h2_o2_ar_burner_flame.csv', quiet=False)
"""
# Get concentration indices
iA4 = gas.species_index('A4')
iC2H2 = gas.species_index('C2H2')
iCO = gas.species_index('CO')
iH = gas.species_index('H')
iH2 = gas.species_index('H2')
iH2O = gas.species_index('H2O')
iO2 = gas.species_index('O2')
iOH = gas.species_index('OH')

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

# Convert pyrene concentration to 10**index
new_concentration = []

for element in f.X[gas.species_index('A4')]:
    new_concentration.append(math.log(element, 10))

# Plot pyrene concentration
conc_comp.plot(c_distance_vals, concentration_vals, 'o', color='g')
conc_comp.plot(distance_grid, fitted_concentration, '-', color='g', label='measured curve')
conc_comp.plot(plot_grid, new_concentration, '-', color='b', label='calculated curve')
conc_comp.title.set_text('Concentration (mole fraction)')
conc_comp.legend()
plt.show()

"""

# Get pyrene concentration and temp
max_index, max_pyrene_molar_concentration = max(enumerate(f.concentrations[iA4, :]), key=operator.itemgetter(1)) #kmoles/m^3
min_index, min_pyrene_molar_concentration = min(enumerate(f.concentrations[iA4,:]), key=operator.itemgetter(1))
pyrene_number_concentration = max_pyrene_molar_concentration * 1000 *  6.0221409e23
pyrene_mole_fraction = f.X[pyrene_index, index]
reaction_temp = f.T[index]

# plot major species
plt.plot(f.grid, f.X[gas.species_index('C4H2'), :], 'b', label='C4H2')
plt.plot(f.grid, f.X[gas.species_index('C6H2'), :], 'g', label='C6H2')
plt.plot(f.grid, f.X[gas.species_index('A1'), :], 'r', label='benzene')
plt.plot(f.grid, f.X[gas.species_index('C2H4'), :], 'm', label='C2H4')
plt.plot(f.grid, f.X[gas.species_index('A2'), :], 'y', label='naphthalene')
plt.xticks(np.arange(0, .1, .02))
plt.yscale('log')
plt.ylim(10e-8, 10e-2)
plt.xlabel('distance (m)')
plt.ylabel('concentration (molefraction)')
plt.title('Calculated mole fractions of major species')
plt.legend(loc='upper right')
plt.show()


# plot minor species
plt.plot(f.grid, f.X[gas.species_index('CO'), :], 'b', label='CO')
plt.plot(f.grid, f.X[gas.species_index('H2'), :], 'g', label='H2')
plt.plot(f.grid, f.X[gas.species_index('C2H2'), :], 'r', label='C2H2')
plt.plot(f.grid, f.X[gas.species_index('H2O'), :], 'm', label='H2O')
plt.plot(f.grid, f.X[gas.species_index('CO2'), :], 'y', label='CO2')
plt.plot(f.grid, f.X[gas.species_index('O2'), :], 'c', label='O2')
plt.plot(f.grid, f.X[gas.species_index('CH4'), :], 'k', label='CH4')
plt.xticks(np.arange(0, .1, .02))
plt.yscale('log')
plt.ylim(10e-4, 0)
plt.xlabel('distance (m)')
plt.ylabel('concentration (molefraction)')
plt.title('Calculated mole fractions of minor species')
plt.legend(loc='upper right')
plt.show()
"""
