# This is an exercise to practice Euler's Method to differential equations
# Credit for exercise is attributed to: https://github.com/afeiguin/comp-phys/blob/master/01_01_euler.ipynb

import numpy as np
import matplotlib.pyplot as plt

# Starting Values
T0 = 10.  # degrees C
Ts = 83.  # degrees C
r = 0.1  # min^-1
t = 0  # sec
dt = 0.05  # time step
tmax = 60.  # max time
nsteps = int(tmax/dt)  # number of steps

# Create numpy arrays to store x and y vals
my_time = np.zeros(nsteps)
my_temp = np.zeros(nsteps)

# Calculate T at each timestep
T = T0
my_temp[0] = T0
for i in range(1,nsteps):
    T = T - r*(T-Ts)*dt
    my_time[i] = i*dt
    my_temp[i]  = T

plt.plot(my_time, my_temp, ls='-', lw=3)
plt.xlabel('time (s)')
plt.ylabel('Temp (C)')
plt.show()

# Challenge: print temp at t=10 as function of delta t
t = 10
T = (r*Ts*t+T0)/(1 + r*t)
dt = (T - T0)/(-1*r*(T-Ts))
print(dt)
