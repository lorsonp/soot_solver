import numpy as np
import math
import random
from matplotlib import pyplot as plt
from IPython.display import clear_output

#  Exercise 10.1
#  Item 1: Find estimate I(N) for the integral of f(x) = 4*np.sqrt(1-x**2) in the interval (0,1)

num_samples=10000


# Define a function to get a random number from a uniform distribution between two input values
def get_random_number(min_val, max_val):
    range = max_val - min_val
    choice = random.uniform(0, 1)
    return min_val + range*choice


# Define integrand function
def f_of_x(x):
    return 4*np.sqrt(1-x**2)

def crude_monte_carlo(num_samples):
    lwr_bnd = 0
    upr_bnd = 1

    sum_of_samples=0
    for i in range(num_samples):
        x = get_random_number(lwr_bnd, upr_bnd)
        sum_of_samples += f_of_x(x)

    return (upr_bnd - lwr_bnd) * float(sum_of_samples/num_samples)


print("The crude Monte Carlo approximation is: %s" % crude_monte_carlo(10000))

# Determine the Variance of the estimation (how much f(x) varies in the domain of x)
def get_crude_MC_variance(num_samples=10000):
    """
    This function returns the variance for the Crude Monte Carlo
    """
    int_max = 1  # the max of our integration range

    # find the average of squares
    running_total = 0
    for i in range(num_samples):
        x = get_random_number(0, int_max)
        running_total += f_of_x(x)**2
    sum_of_sqs = (int_max*running_total/num_samples)

    # find square of ave
    running_total = 0
    for i in range(num_samples):
        x = get_random_number(0, int_max)
        running_total = f_of_x(x)
    sq_ave = (int_max*running_total/num_samples)**2

    return sum_of_sqs - sq_ave

crude_MC_variance = get_crude_MC_variance()
crude_MC_SE = np.sqrt(crude_MC_variance/num_samples)

print("Variance: %s" % crude_MC_variance)
print("Standard Error: %s" % crude_MC_SE)
