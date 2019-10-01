# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script defines variables for the modified species equation             #                                                                            #
#                                                                             #
# Content for this code was derived from a Soot subroutine originally         #
# developed by Kenneth Revzan, Nancy Brown, and Michael Frenklach.            #
#                                                                             #
# Credit to the source code is given to:                                      #
# http://combustion.berkeley.edu/soot/codes/codes.html.                       #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
import math

# Primary Equation
# (d/dz) (-Di_0 * (dMr_2_3/dz) + vT*Mr + rho*v * (d(Mr/rho)/dz) - Qr = 0, r = 1, ..., Nmoments, Nmoments = 6

# Constants
rho = 1.8e3                        # density of soot (kg/m^3)
Avo = 6.022e23                     # Avogadro's number (particles/mole)
k_B = 1.38e-23                     # Boltzmann's constant (m^2kg/(s^2K))
T_max = 1992                       # Max flame temperature (K)
m_C = 2.004e-23                    # Mass of 1 carbon atom (kg)
velocity = 0.204                   # m/s
mean_mass = (0.55*39.95 + 0.2143*32 + 0.2357*26.04)/(1000*Avo)  # kg
# Note: mean mass of gas has been calculated for initial mole fractions, whereas it may need to be calculated dynamically in the flame simulation

# Di_0 : Diffusion rate of smallest soot particle (K*0.5/s)
Di_0 = 3/(2*rho) * (mean_mass*k_B*T/(2*np.pi))**(1/2) * (1 + np.pi*alpha/8)**(-1) * 1/di_0**2

# T : temperature. Cantera will provide this info. (K)

# alpha : thermal accommodation coefficient (no units)
a = 12.65 - 0.00563*T
b = -1.38 + 0.00068*T
mu_1 = 32                          # 1st size moment = average particle size (# of carbon atoms)
alpha = np.tanh(a/math.log10(mu_1) + b)
# Note: For additional precision, alpha is computed as a function of temperature.
# If you wish to simplify this, you can use a constant instead: alpha = 0.587

# di_0 : Diameter of smallest soot particle (m)
di_0 = (6*m_C/(np.pi*rho))**(1/3) * (32)**(1/3)

# v_T : Thermal diffusion velocity (m/s)
v_T = -3/4 * (1 + np.pi*alpha/8)**-1 * gas_viscosity/(rho*T) * (dT/dz)
# Note: Not sure how to adjust this equation to account for the (dT/dz)
