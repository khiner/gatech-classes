import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt

# Constants
NTl = 1000            # Number of thallium atoms
NPb = 0               # Number of lead atoms
tau = 183.2           # Half life of thallium in seconds
h = 1.0               # Size of time-step in seconds
p = 1 - 2**(-h/tau)   # Probability of decay in one step
tmax = 1000           # Total time

# Lists of plot points
tpoints = np.arange(0.0,tmax,h)
Tlpoints = []
Pbpoints = []

# Main loop
for t in tpoints:
    Tlpoints.append(NTl)
    Pbpoints.append(NPb)

    # Complete in-class
    # Calculate the number of atoms that decay in this timestep
    #

# Make the graph
plt.plot(tpoints,Tlpoints)
plt.plot(tpoints,Pbpoints)
plt.xlabel("Time")
plt.ylabel("Number of atoms")
plt.show()
