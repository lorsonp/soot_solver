# This is a Monte Carlo exercise to practice the crude and importance sampling MC methods
# Credit for exercise is attributed to: https://www.ias.ac.in/article/fulltext/reso/019/08/0713-0739

import numpy as np
import math
import random
import matplotlib.pyplot as plt

PI = 3.1415926
a = 0
b = PI
num_samples = 10000

# Integrand
def f_of_x(x):
    return 1.0/(x**2 + (math.cos(x))**2)


# Random number generator
def get_rand_num(min_val, max_val):
    return random.uniform(min_val, max_val)


# Crude Monte Carlo Evaluation
def crude_MC(a, b, num_samples):

    sum_of_samples = 0
    for i in range(num_samples):
        x = get_rand_num(a, b)
        sum_of_samples += f_of_x(x)

    return float(sum_of_samples * (b - a)/num_samples)


print("The crude Monte Carlo approximation is: " + str(crude_MC(0, PI, 10000)))


# Crude MC variance
def crude_MC_var(a, b, num_samples):

    # find ave of squares
    running_total = 0
    for i in range(num_samples):
        x = get_rand_num(a, b)
        running_total += f_of_x(x)**2
    ave_of_squares = running_total * (b - a)/num_samples


    # find square of ave
    running_total = 0
    for i in range(num_samples):
        x = get_rand_num(a, b)
        running_total += f_of_x(x)
    square_of_ave = (running_total * (b - a)/num_samples)**2

    return math.fabs(ave_of_squares - square_of_ave)

print("The crude Monte Carlo variance is: " + str(crude_MC_var(0, PI, 10000)))


# Importance function
def g_of_x(x, A, lamda):
    return A * math.exp(-1 * lamda * x)


# Inverse function
def G_inv_of_r(lamda):
    return (-1/lamda) * math.log((1 - math.exp(-PI * lamda))/lamda)


# Importance Sampling MC Variance
def is_MC_var(a, b, lamda, num_samples):
    A = lamda / (1 - math.exp(-1*PI*lamda))

    # find ave of squares
    running_total = 0
    for i in range(num_samples):
        x = get_rand_num(a, b)
        running_total += (f_of_x(x)/g_of_x(A, lamda, x)) ** 2

    ave_of_squares = running_total / num_samples

    # find square of ave
    running_total = 0
    for i in range(num_samples):
        x = get_rand_num(a, b)
        running_total += f_of_x(x)/g_of_x(A, lamda, x)
    square_of_ave = (running_total / num_samples) ** 2

    return ave_of_squares - square_of_ave


# Find lowest variance by testing a range of lambda values
test_lamdas = np.arange(0.05, 1.60, 0.05)
variance = np.zeros(len(test_lamdas))

for i in range(len(test_lamdas)):
    variance[i] = is_MC_var(a=0, b=PI, lamda=test_lamdas[i], num_samples=10000)

plt.figure()
plt.plot(test_lamdas, variance, 'o')
plt.show()

min_var = np.where(variance == variance.min())
print("Using importance sampling, the minimum variance is found to be " + str(variance[16]) + " when lambda is " + str(test_lamdas[16]) + ".")


# Run importance sampling MC using optimal lambda value
def is_MC(a=0, b=PI, lamda=test_lamdas[16], num_samples=10000):
    A = lamda / (1 - math.exp(-1 * PI * lamda))

    running_total = 0
    for i in range(num_samples):
        r = get_rand_num(0, 1)
        running_total += f_of_x(G_inv_of_r(lamda=test_lamdas[16])) / g_of_x(G_inv_of_r(lamda=test_lamdas[16]), A, lamda=test_lamdas[16])

    return float(running_total/num_samples)


print("The importance sampling Monte Carlo approximation is: " + str(is_MC()) + ".")

