%matplotlib inline
from random import random,randrange
from math import exp,pi
import numpy as np
import matplotlib.pyplot as plt

T = 10.0
N = 1000
steps = 250000

# Create a 2D array to store the quantum numbers
n = np.ones([N,3],int)

# Main loop
eplot = []
E = 3*N*pi*pi/2
for k in range(steps):

    # Choose the particle and the move
    i = randrange(N)
    j = randrange(3)
    
    # Complete in-class
    #
    # Determine a random move set and associate change in energy
    #

    # Complete in-class
    #
    # Decide whether to accept the move. If so, change the state and total energy
    #

    eplot.append(E)

# Make some graphs
plt.plot(eplot)
plt.xlabel("Iteration")
plt.ylabel("Energy")
plt.show()

plt.hist(n, 3, label=['x','y','z'])
plt.legend(loc='best')
plt.show()