# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This script includes modifications to be made to Cantera.                   #
#                                                                             #
# The modifications will allow Cantera to return the "rates" of the first six #
# moments of the soot PSDF along with adjustments to the production rate of   #
# pyrene to account for pyrene consumption during inception.                  #
#                                                                             #
# This code was modified from a Soot subroutine originally developed by       #
# Kenneth Revzan, Nancy Brown, and Michael Frenklach.                         #
#                                                                             #
# Credit to the source code is given to:                                      #
# http://combustion.berkeley.edu/soot/codes/codes.html.                       #
#                                                                             #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# --------------- #
# required inputs #
# --------------- #

PAH_concentration                           # concentration of PAH (mol/cm3)
Y = [PAH_concentration, 0, 0, 0, 0, 0, 0]   # array of chemical species concentrations (mol/cm3) and soot moments (cm3)
T                                           # temperature (K)
P                                           # pressure (atm)
N_moments = 6                               # number of soot moments to be calculated
iPAH = 0                                    # index of pyrene in array Y
iM0 = 2                                     # index of the first moment in array Y
freeMolOnly = True                          # coagulation rates are calculated for free-molecular regime only
M0_cutoff = 0                               # below this value, only nucleation rates are calculated
CPerPAH = 16                                # number of carbon atoms in pyrene
diamPAH = 7.90369e-8                       # diameter of pyrene (cm)

# ------------------------------------------------ #
# Variables/Equations needed for rate calculations #
# ------------------------------------------------ #

CPerDimer = 2 * CPerPAH
AVOGADRO = 6.024e23
PI = 3.14159
BOLTZMANN = 1.3807e-16                      # m^2*g/(s^2*K)
C_MASS = 12 * 1.67e-24                      # g
RHO_SOOT = 1.8                              # g/cm^3
cbNucl = 2.2 * (4 * PI * BOLTZMANN / (
            C_MASS * CPerPAH)) ** 0.5 * diamPAH ** 2 * AVOGADRO ** 2    # pyrene coagulation kernel
CB0 = 2.2 * (3 * C_MASS / (4 * PI * RHO_SOOT)) ** (1 / 6) * (
            6 * BOLTZMANN / RHO_SOOT) ** 0.5                            # soot coagulation kernel

# --------------------------- #
# Initialize output variables #
# --------------------------- #

nucRate = []
coagRate = []

for index in range(0, 6):
    nucRate(index) = 0
    coagRate(index) = 0

prodPAH = 0

# --------------------------- #
# Nucleation rate for moments #
# --------------------------- #

if Y(iPAH) > 0:
    nucRate(0) = cbNucl * T ** 0.5 * Y(iPAH) * Y(iPAH)
else:
    nucRate(0) = 0

for index in range(1, 7):
    nucRate(index) = CPerDimer * nucRate(index - 1)

# --------------------------------------- #
# PAH consumption rate due to nucleation) #
# --------------------------------------- #

prodPAH = -2 * nucRate(0) / AVOGADRO

# ---------------------------- #
# Coagulation rate for moments #
# ---------------------------- #

coagRate = coagulation(P, T, M0=Y(iM0), coagRate=0)

# -------------------- #
# Rates of the moments #
# -------------------- #

for index in range(0, 6):
    rate(index) = nucRate(index) + coagRate(index)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# NOTE: The following functions are in FORTRAN syntax.                  #
# They will need to be modified in order to be compatible with Cantera. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# --------------------------------------------------------------------- #
"""
subroutine coagulation(P, T, M0, rate)

# Coagulation rates calculated using Method II of M.Frenklach and S.Harris, J.Colloid Interface Sci., v .118, n .1, 1987,
# pp. 252 - 261. See also M.Frenklach and H.Wang, "Soot Formation in Combustion: Mechanisms and Models" (H. Bockhorn, Ed.),
# Springer - Verlag, Heidelberg, 1994, p. 165.

# Free molecular regime

    CB0_sT = CB0 * T ** 0.5
    
    do n = 1, 4
        f(n) = gridFun(n - 1, 0, 0)
    enddo
    crk = CB0_sT * fitpos(0.5, 4, f)
    do n = 1, 4
        f(n) = gridFun(n - 1, 1, 1)
    enddo
    crk2 = CB0_sT * fitpos(0.5, 4, f)
    do n = 1, 4
        f(n) = gridFun(n - 1, 1, 2)
    enddo
    crk3 = CB0_sT * fitpos(0.5, 4, f)
    do n = 1, 3
        f(n) = gridFun(n - 1, 1, 3)
    enddo
    crk4a = CB0_sT * fitpos(0.5, 3, f)
    do n = 1, 4
        f(n) = gridFun(n - 1, 2, 2)
    enddo
    crk4b = CB0_sT * fitpos(0.5, 4, f)
    do n = 1, 2
        f(n) = gridFun(n - 1, 1, 4)
    enddo
    crk5a = CB0_sT * fitpos(0.5, 2, f)
    do n = 1, 3
        f(n) = gridFun(n - 1, 2, 3)
    enddo
    crk5b = CB0_sT * fitpos(0.5, 3, f)
    
    M02 = M0 * M0
    rate(0) = -0.5 * crk * M02
    rate(1) = 0.
    rate(2) = crk2 * M02
    rate(3) = 3 * crk3 * M02
    rate(4) = (4 * crk4A + 3 * crk4B) * M02
    rate(5) = (5 * crk5A + 10 * crk5B) * M02

return rate

# --------------------------------------------------------------------- #

function gridFun(k, n, m)

# This function used in calculation of free - molecular collision rates, which, for any(n, m), are determined by Lagrange
# interpolation using two or more gridFun(k, n, m)

    implicit none
    
    include 'soot.h'
    
    integer k, n, m
    
    doubleprecision mu(LO_FRAC: HI_FRAC)
    common / MOMENTS / mu
    
    integer i, j, l
    
    gridFun = 0.0
    dol = 0, k
        i = 6 * (k - l + n)
        j = 6 * (l + m)
        gridFun = gridFun + binom(l, k) * (
            & mu(i+1) * mu(j - 3) +
            & 2 * mu(i - 1) * mu(j - 1) +
            & mu(i - 3) * mu(j + 1) )
    enddo

return

end

# --------------------------------------------------------------------- #

function fitpos(point, numTerms, y)

# Perform Lagrange interpolation at `point`, using values in y
# x is in common

    implicit none
    
    include 'soot.h'
    
    doubleprecision point
    integer numTerms
    doubleprecision y(numTerms)
    
    doubleprecision a, prod
    integer n
    
    prod = 1.0
    do n = 1, numTerms
        prod = prod * (point - x(n))
    enddo
    
    fitpos = 1.0
    do n = 1, numTerms
        a = prod / ((point - x(n)) * prime(numTerms, n))
        fitpos = fitpos * (y(n) ** a)
    enddo

return

end
"""
