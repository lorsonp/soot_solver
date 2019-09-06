# This is an exercise which calculates chi squared for exoplanet transits
# Credit for the code is attributed to: https://star.pst.qub.ac.uk/wiki/doku.php/users/kpoppenhaeger/phy1024/start

import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0, 100)
transit_expected = np.ones(len(t))
transit_expected[40:60] = 0.9

measurements = np.ones(len(t))
measurements[47:67] = 0.9
np.random.seed(111)
measurements = measurements + np.random.normal(0, 0.015, len(t))
errors = 0.01 + np.random.uniform(0, 0.01, len(t))

plt.figure(figsize=(8,6))
plt.plot(t, transit_expected)
plt.errorbar(t, measurements, errors, fmt='r.')
plt.show()

chi_sq = ((measurements - transit_expected)**2/errors).sum()