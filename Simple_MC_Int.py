# Exercise in calculating MC integration approximation using the 'Simple' Method: I(N) = (b-a)*<f>
# Credit for source code: https://github.com/afeiguin/comp-phys/blob/master/10_01_montecarlo_integration.ipynb

import numpy as np
import matplotlib.pyplot as plt

# Simple Monte Carlo Integration
ngroups = 16

I = np.zeros(ngroups)
N = np.zeros(ngroups)
E = np.zeros(ngroups)

n0 = 100
for i in range(ngroups):

    N[i] = n0
    r = np.random.random(n0)
    I[i] = 0.
    for j in range(n0):
        x = r[j]
        I[i] += np.sqrt(1 - x ** 2)

    I[i] *= 4. / float(n0)
    E[i] = abs(I[i] - np.pi)
    print
    n0, I[i], E[i]
    n0 *= 2

plt.figure()
plt.plot(N, E, ls='-', c='red', lw=3);
plt.plot(N, 0.8 / np.sqrt(N), ls='-', c='blue', lw=3);
plt.xscale('log')
plt.yscale('log')
plt.show()

n0 = 100000
I = np.zeros(n0)
r = np.random.random(n0)
for j in range(n0):
    x = r[j]
    I[j] = 4 * np.sqrt(1 - x ** 2)


def group_measurements(ngroups):
    global I, n0

    nmeasurements = n0 / ngroups
    for n in range(ngroups):
        Ig = 0.
        Ig2 = 0.
        for i in range(int(n * nmeasurements), int((n + 1) * nmeasurements)):
            Ig += I[i]
            Ig2 += I[i] ** 2
        Ig /= nmeasurements
        Ig2 /= nmeasurements
        sigma = Ig2 - Ig ** 2
        print(Ig, Ig2, sigma)


group_measurements(10)
print("=============================")
group_measurements(20)
print("=============================")
group_measurements(1)