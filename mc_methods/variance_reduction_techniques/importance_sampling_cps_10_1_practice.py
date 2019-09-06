# Exercise in MC integration approximation using importance sampling
# Credit for source code: https://github.com/afeiguin/comp-phys/blob/master/10_01_montecarlo_integration.ipynb

import numpy as np
import matplotlib.pyplot as plt
import math
import random

# plot functions in question
plt.xlim(0, 10)
x = np.arange(0, 10, 0.1)
plt.plot(x, np.exp(-1*x))
plt.plot(x, x**1.5*np.exp(-x))
plt.show()

# Trapezoidal Integration
def trapezoids(func, xmin, xmax, nmax):
    Isim = func(xmin) + func(xmax)
    h = (xmax - xmin)/nmax
    for i in range(1, nmax):
        x = xmin+i*h
        Isim += 2*func(x)

    Isim *= h/2
    return Isim


# function
def f_of_x(x):
    return x**1.5*np.exp(-x)


print("Trapezoids: ", trapezoids(f_of_x, 0., 20., 100000))

# Simple MC
n0 = 100000
r = np.random.random(n0)
r = 20*r

Itot = np.sum((r)**1.5*np.exp(-r))
print("Simple Monte Carlo: ", Itot*20/n0)


# Importance Sampling
r = np.random.random(n0)

x = -np.log(r)
Itot = np.sum(x**1.5)
print("Importance Sampling: ", Itot/n0)

# plot next practice function
fig, ax = plt.subplots()
plt.xlim(0, np.pi)
x = np.arange(0, np.pi, 0.05)
ax.plot(x, 1./(x**2 + np.cos(x)**2), label='Integrand')
ax.plot(x, np.exp(-x), label='Weight Function 1')
ax.plot(x, np.exp(-2*x), label='Weight Function 2')
plt.title('Importance Sampling for Reduced Variance')
legend = ax.legend(loc=0)
plt.show()

# function
def g_of_x(x):
    return 1./(x**2 + np.cos(x)**2)


print("Trapezoids: ", trapezoids(g_of_x, 0., np.pi, 100000))

# Simple MC
a = np.arange(0.1, 2.1, 0.1)
I = np.arange(0.1, 2.1, 0.1)

r = np.random.random(n0)
r = np.pi*r

I0 = np.sum(1./(r**2+np.cos(r)**2))
print("Simple Monte Carlo: ", I0*np.pi/n0)

# Importance Sampling
print("Importance Sampling: ")
x = -np.log(r)
