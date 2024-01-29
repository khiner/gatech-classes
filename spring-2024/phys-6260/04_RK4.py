import numpy as np
import matplotlib.pyplot as plt

def f(r,t):
    x = r[0]
    y = r[1]
    fx = x*y - x
    fy = y - x*y + np.sin(t)**2
    return np.array([fx,fy],float)

a = 0.0
b = 10.0
N = 1000
h = (b-a)/N

tpoints = np.arange(a,b,h)
r = np.zeros((N,2))
r[0] = [1.0, 1.0]

# To be completed in class: loop over time to solve for x and y, which
# will be stored in r[]


plt.plot(tpoints, r[:,0], label='x')
plt.plot(tpoints, r[:,1], label='y')
plt.xlabel("t")
plt.ylabel("x,y")
plt.show()
