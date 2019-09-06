# This is an exercise to practice methods of integration as applied to laws of motion
# Credit for exercise is attributed to: https://github.com/afeiguin/comp-phys/blob/master/01_03_eqs_of_motion.ipynb

import numpy as np
import matplotlib.pyplot as plt


# Challenge 1.2: Write a program to solve the 1d equations of motion for a falling object.

# define a class "particle"
class particle(object):

    def __init__(self, mass=1., y=0, v=0):
        self.mass = mass
        self.y = y
        self.v = v

    def euler(self, f, dt):
        self.y = self.y + self.v*dt
        self.v = self.v + f/self.mass*dt

    def euler_cromer(self, f, dt):
        self.v = self.v + f/self.mass*dt
        self.y = self.y + self.v*dt

# input variables
g = 9.8
mass = 0.01
y0 = 300.
v0 = 0.
vt = 30.

dt = 0.5

gforce = g*mass

p = particle(mass, y0, v0)

y = [y0]
v = [v0]
t = [0.]

while p.y > 0.:
    fy = -gforce
    p.euler(fy, dt)
    y.append(p.y)
    v.append(p.v)
    t.append(t[-1]+dt)

t_data = np.array(t)  # convert list to array for plotting

