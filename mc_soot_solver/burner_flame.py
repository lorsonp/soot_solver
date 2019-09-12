"""
A burner-stabilized lean premixed hydrogen-oxygen flame at low pressure.
"""

import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import operator

p = 0.12 * ct.one_atm
tburner = 300.0
reactants = 'AR:0.55, O2:0.2143, C2H2:0.2357'  # premixed gas composition
width = 0.1 # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('abf_mech90torr.cti')
gas.TPX = tburner, p, reactants

mdot = 0.204 * gas.density  # kg/(m^2*s)
print(gas.density)
print(mdot)

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show_solution()

f.transport_model = 'Mix'
f.solve(loglevel, auto=True)
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

"""
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


# plot pyrene
plt.plot(f.grid, f.X[iA4, :], 'b', label='pyrene')
plt.xticks(np.arange(0, .05, .01))
plt.yscale('log')
plt.ylim(10e-9, 10e-5)
plt.xlabel('distance (m)')
plt.ylabel('concentration (molefraction)')
plt.title('Calculated concentration of pyrene')
plt.show()
"""

