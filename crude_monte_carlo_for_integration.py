# This is a Monte Carlo exercise which estimates the value of an integral by multiplying the average value of the integrand by the domain of integration.
# Variance and standard error are then calculated.
# Credit for the code is attributed to: https://towardsdatascience.com/monte-carlo-simulations-with-python-part-1-f5627b7d60b0

import numpy as np
import random

PI = 3.1415926
e = 2.71828
num_samples = 10000


# Define a function to get a random number from a uniform distribution between two input values
def get_random_number(min_val, max_val):
    range = max_val - min_val
    choice = random.uniform(0, 1)
    return min_val + range*choice


# Define integrand function
def f_of_x(x):
    # x must be in radians
    return (e**(-1*x))/(1 + (x - 1)**2)


def crude_monte_carlo(num_samples):
    """
    This function performs the crude Monte Carlo for f(x) from (0,5)
    """
    lower_bnd = 0
    upper_bnd = 5

    sum_of_samples = 0
    for i in range(num_samples):
        x = get_random_number(lower_bnd, upper_bnd)
        sum_of_samples += f_of_x(x)

    return (upper_bnd - lower_bnd) * float(sum_of_samples/num_samples)


# Determine the Variance of the estimation (how much f(x) varies in the domain of x)
def get_crude_MC_variance(num_samples=10000):
    """
    This function returns the variance for the Crude Monte Carlo
    """
    int_max = 5  # the max of our integration range

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
