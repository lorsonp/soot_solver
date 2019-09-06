#  Reimann Sums
#  credit: Larry Feldman
# https://www.youtube.com/watch?v=b8BvG5MKpzU

#  rightsum = y[1](x[1] - x[0]) + . . . + y[n-1](x[n-1] - x[n-2])
#  leftsum = y[0](x[1] - x[0]) + . . . + y[n-2](x[n-1] - x[n-2])
#  trapsum = (rightsum + leftsum)/2

import sys

xvals = (0, 1, 2, 3, 4, 5)
yvals = (0, 1, 4, 9, 16, 25)  # y = x^2


if len(xvals) != len(yvals):
    print("The number of x and y values must be the same")
    sys.exit()


rightsum = 0
leftsum = 0


for i in range(1, len(xvals), 1):
    rightsum = rightsum + yvals[i]*(xvals[i] - xvals[i-1])
    leftsum = leftsum + yvals[i - 1]*(xvals[i] - xvals[i-1])


trapsum = (leftsum + rightsum)/2.
print("Right sum: " + str(rightsum))
print("Leftt sum: " + str(leftsum))
print("Trap sum: " + str(trapsum))
