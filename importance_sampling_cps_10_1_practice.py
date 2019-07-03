# Exercise in MC integration approximation using importance sampling
# Credit for source code: https://github.com/afeiguin/comp-phys/blob/master/10_01_montecarlo_integration.ipynb

import numpy as np
import matplotlib.pyplot as plt
import math
import random


# function
def f_of_x(x):
    return math.exp(-1*x**2)


f2 = np.vectorize(f_of_x)

x = np.arange(0, 1, 0.001)
plt.figure()
plt.plot(x, f2(x))
plt.show()
