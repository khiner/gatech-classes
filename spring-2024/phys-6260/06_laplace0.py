import numpy as np
import matplotlib.pyplot as plt

# Constants
M = 100         # Grid squares on a side
V = 1.0         # Voltage at top wall
target = 3e-4   # Target accuracy

# Create arrays to hold potential values
phi = np.zeros([M+1,M+1],float)
phi[0,:] = V  # top wall
phiprime = np.zeros([M+1,M+1],float)

# Main loop
delta = 1.0
iteration = 0
while delta>target:

    # Calculate new values of the potential
    for i in range(M+1):
        for j in range(M+1):
            #
            # Complete in-class for phiprime, the updated solution
            # Remember to keep the boundaries at zero
            #

    # Calculate maximum difference from old values
    delta = np.max(abs(phi-phiprime))
    if iteration % 10 == 0:
        print("Iteration %3d: max. residual = %12.6g" % (iteration, delta))

    # Swap the two arrays around
    phi,phiprime = phiprime,phi
    iteration += 1

# Make a plot
plt.imshow(phi)
plt.gray()
plt.show()
