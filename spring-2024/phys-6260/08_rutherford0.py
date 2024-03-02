import numpy as np
from numpy.random import random

# Constants
Z = 79
e = 1.602e-19
E = 7.7e6*e
epsilon0 = 8.854e-12
a0 = 5.292e-11
sigma = a0/100
N = 100000

# Complete in-class
# Function to generate two Gaussian random numbers
# Optional: Use array arithmetic to speed-up
def gaussian():
    x = y = 0
    return x,y

count = 0
for i in range(N):
    # Complete in-class
    # For every particle, calculate a random position (x,y) and determine 
    # whether it's back-scattered (b < some critical radii)

print("%d particles were reflected out of %d" % (count,N))
