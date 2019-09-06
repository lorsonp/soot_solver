# This is a Monte Carlo exercise which estimates the value of an integral by multiplying the average value of the integrand by the domain of integration.
# Variance and standard error are calculated and then importance sampling is employed.
# Credit for source code is attributed to: https://towardsdatascience.com/monte-carlo-simulations-with-python-part-1-f5627b7d60b0

import numpy as np
import math
import random
import matplotlib.pyplot as plt
from IPython.display import clear_output

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


print("MC Estimate: " + str(crude_monte_carlo(num_samples)))


# Determine the Variance of the estimation (how much f(x) varies in the domain of x)
def get_crude_mc_variance(num_samples=10000):
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
        running_total += f_of_x(x)
    sq_ave = (int_max*running_total/num_samples)**2

    return math.fabs(sum_of_sqs - sq_ave)


crude_mc_variance = get_crude_mc_variance()
crude_mc_SE = np.sqrt(crude_mc_variance/num_samples)

print("Variance: %s" % crude_mc_variance)
print("Standard Error: %s" % crude_mc_SE)


# this is the template of our weight function g(x)
def g_of_x(x, A, lamda):
    e = 2.71828
    return A * math.pow(e, -1 * lamda * x)


def inverse_G_of_r(r, lamda):
    return (-1 * math.log(float(r))) / lamda


def get_IS_variance(lamda, num_samples):
    """
    This function calculates the variance if a Monte Carlo
    using importance sampling.
    Args:
    - lamda (float) : lamdba value of g(x) being tested
    Return:
    - Variance
    """
    A = lamda
    int_max = 5

    # get sum of squares
    running_total = 0
    for i in range(num_samples):
        x = get_random_number(0, int_max)
        running_total += (f_of_x(x) / g_of_x(x, A, lamda)) ** 2

    sum_of_sqs = running_total / num_samples

    # get squared average
    running_total = 0
    for i in range(num_samples):
        x = get_random_number(0, int_max)
        running_total += f_of_x(x) / g_of_x(x, A, lamda)
    sq_ave = (running_total / num_samples) ** 2

    return sum_of_sqs - sq_ave


# get variance as a function of lambda by testing many
# different lambdas

test_lamdas = [i * 0.05 for i in range(1, 61)]
variances = []

for i, lamda in enumerate(test_lamdas):
    print(f"lambda {i + 1}/{len(test_lamdas)}: {lamda}")
    A = lamda
    variances.append(get_IS_variance(lamda, 10000))
    clear_output(wait=True)

optimal_lamda = test_lamdas[np.argmin(np.asarray(variances))]
IS_variance = variances[np.argmin(np.asarray(variances))]

print(f"Optimal Lambda: {optimal_lamda}")
print(f"Optimal Variance: {IS_variance}")
print(f"Error: {(IS_variance / 10000) ** 0.5}")

# Run Simulation
def importance_sampling_MC(lamda, num_samples):
    A = lamda

    running_total = 0
    for i in range(num_samples):
        r = get_random_number(0, 1)
        running_total += f_of_x(inverse_G_of_r(r, lamda=lamda)) / g_of_x(inverse_G_of_r(r, lamda=lamda), A, lamda)
    approximation = float(running_total / num_samples)
    return approximation


# run simulation
num_samples = 10000
approx = importance_sampling_MC(optimal_lamda, num_samples)
variance = get_IS_variance(optimal_lamda, num_samples)
error = (variance / num_samples) ** 0.5

# display results
print(f"Importance Sampling Approximation: {approx}")
print(f"Variance: {variance}")
print(f"Error: {error}")
