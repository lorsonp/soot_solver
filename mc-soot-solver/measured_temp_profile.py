import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import plotly

distance = [0, 3.81, 6.30, 12.86, 20.53, 29.21, 38.52, 48.15, 58.89, 70.74, 82.96]  # mm
measured_temp = [300, 1945.52, 1992, 1950, 1900, 1850, 1800, 1750, 1700, 1650, 1600]  # K
# splines = interpolate.splrep(distance, measured_temp)

range_1 = round(1000*(3.81/90))
range_2 = round(1000*(6.3-3.819)/90)
range_3 = round(1000*((90-6.309)/90))

x_vals_1 = np.linspace(0, 3.81, range_1)
x_vals_2 = np.linspace(3.8, 6.3, range_2)
x_vals_3 = np.linspace(6.309, 90, range_3)


# spline_y_vals = interpolate.splev(x_vals_1, splines)
x_1 = np.array([0, 3.81])
y_1 = np.array([300, 1945.52])

x_2 = np.array([3.81, 6.30])
y_2 = np.array([1945.52, 1992])

f_1 = np.poly1d(np.polyfit(x_1, y_1, 1))
f_2 = np.poly1d(np.polyfit(x_2, y_2, 2))
f_3 = np.poly1d(np.polyfit(distance[2:-1], measured_temp[2:-1], 2))

plt.plot(distance, measured_temp, 'o', color='r')
plt.plot(x_vals_1, f_1(x_vals_1), '-', color='g')
plt.plot(x_vals_2, f_2(x_vals_2), '-', color='g')
plt.plot(x_vals_1, f_1(x_vals_1), '-', color='g')
plt.plot(x_vals_3, f_3(x_vals_3), '-', color='g')
plt.show()
