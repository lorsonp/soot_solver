# This is an exercise which calculates goodness of fit for exoplanet transits
# Credit for the code is attributed to: https://star.pst.qub.ac.uk/wiki/doku.php/users/kpoppenhaeger/phy1024/start

import numpy as np
import matplotlib.pyplot as plt

t = np.arange(0, 30)
# expected transit model
transit_expected = np.ones(len(t))
# setting expected bottom of transit
transit_expected[5:10] = 0.9

y = np.ones(len(t))
# setting bottom of transit
y[15:20] = 0.9
# add random noise
np.random.seed(123456)
y = y + np.random.normal(0, 0.005, len(t))
# define error bars for data
e = 0.005 + np.random.normal(0, 0.001, len(t))

plt.figure()
plt.plot(t, transit_expected, '-')
plt.errorbar(t, y, e, fmt = 'r.')
plt.axis([0, 30, 0.8, 1.10])
plt.show()