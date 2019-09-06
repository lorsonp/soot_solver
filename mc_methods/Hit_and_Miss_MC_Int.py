# Exercise in calculating MC integration approximation using the 'Hit and Miss' Method: I(N) = (Nin/N)*H*(b-a)
# Credit for source code: https://github.com/afeiguin/comp-phys/blob/master/10_01_montecarlo_integration.ipynb

import numpy as np
import matplotlib.pyplot as plt

# Hit and miss Monte Carlo integration
ngroups = 16

I = np.zeros(ngroups)
N = np.zeros(ngroups)
E = np.zeros(ngroups)

n0 = 100
for i in range(ngroups):

    N[i] = n0
    x = np.random.random(n0)
    y = np.random.random(n0)
    I[i] = 0.
    Nin = 0
    for j in range(n0):
        if (y[j] < np.sqrt(1 - x[j] ** 2)):
            Nin += 1

    I[i] = 4. * float(Nin) / float(n0)
    E[i] = abs(I[i] - np.pi)
    print(n0, Nin, I[i], E[i])
    n0 *= 2

plt.figure()
plt.plot(N, E, ls='-', c='red', lw=3);
plt.plot(N, 0.8 / np.sqrt(N), ls='-', c='blue', lw=3);

plt.xscale('log')
plt.yscale('log')
plt.show()
