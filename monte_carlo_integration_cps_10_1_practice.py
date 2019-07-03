import numpy as np
import math
import random
from matplotlib import pyplot as plt
from IPython.display import clear_output


# Define a function to get a random number from a uniform distribution between two input values
def get_random_number(min_val, max_val):
    range = max_val - min_val
    choice = random.uniform(0, 1)
    return min_val + range*choice


#  Exercise 10.1
#  Item 1: Find estimate I(N) for the integral of f(x) = 4*np.sqrt(1-x**2) in the interval (0,1)

a = 0
b = 1
PI = 3.1415926

# Define integrand function
def f_of_x(x):
    return 4*np.sqrt(1-x**2)


def crude_monte_carlo(a, b, num_samples):

    sum_of_samples=0
    for i in range(num_samples):
        x = get_random_number(a, b)
        sum_of_samples += f_of_x(x)

    return float((b-a)*sum_of_samples/num_samples)


print("The crude Monte Carlo approximation for 10,000 samples is: %s" % crude_monte_carlo(a, b, num_samples=10000))


# Determine difference between approximation and exact value of PI (10,000 samples)
Error_10k = abs(PI - crude_monte_carlo(a, b, num_samples=10000))
print("The error for 10,000 samples is: %s" % Error_10k)

# Make loglog plot (error as function of N)

x = np.arange(1000, 100000, 1000)
error = np.zeros(len(x))

for i in range(len(x)):
    error[i] = abs(PI - crude_monte_carlo(a=0, b=1, num_samples=x[i]))

plt.figure()
plt.loglog(x, error)
plt.show()

