# This is a Monte Carlo exercise which simulates a probability problem: Given two six-sided die, what is the chance of rolling 11 or higher?
# A series of simulations are then plotted and compared to the true value, and standard deviation is calculated.
# Credit for the code is attributed to: https://star.pst.qub.ac.uk/wiki/doku.php/users/kpoppenhaeger/phy1024/start

import numpy as np
import matplotlib.pyplot as plt
import statistics as st


def montecarlodice(myseed):  # Dice Simulation Function
    length = 10000  # how many random numbers per single die
    np.random.seed(myseed)  # set seed for the first die
    die_1 = np.random.random_integers(1, 6, length)
    np.random.seed(myseed + 1)  # set different seed for second die
    die_2 = np.random.random_integers(1, 6, length)
    result = die_1 + die_2
    mymask = result >= 11
    montecarloresult = np.float(mymask.sum())/np.float(len(result))
    return montecarloresult


myseedlist = np.array([100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])  # array of 10 different seeds
myresults = np.zeros(10)  # empty array to store results of 10 simulations

# make a loop to run 10 simulations
for i in np.arange(0, 10):
    myresults[i] = montecarlodice(myseedlist[i])

print(myresults)

# plot the values and see how they compare to the true value
true_result = 1./12.

plt.figure()
x = np.arange(0, len(myresults))
plt.plot(x, myresults, 'o')  # plot simulation results as dots
plt.plot(x, np.ones(len(x))*true_result, '-')  # plot the true result as a line
plt.show()

print("Standard Deviation is: " + str(st.stdev(myresults, true_result)))
