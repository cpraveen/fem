import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

# Plot lagrange polynomials
for n in range(7,1,-1):
    plt.figure()
    xs = np.linspace(0.0, 1.0, n)
    for i in range(n):
        values = np.zeros((n,1))
        values[i] = 1.0
        fun = lagrange(xs, values)
        x = np.linspace(0,1,100)
        plt.plot(x, fun(x))
        plt.plot(xs,0*xs,'o')
        plt.xlabel('x')
        plt.grid(True)
        plt.title('Number of points = '+str(n))

# Plot the points
plt.figure()
for n in range(2, 8):
    x = np.linspace(0.0, 1.0, n)
    plt.plot(x, (10-n)*np.ones(n), 'o-')
plt.yticks([])
plt.title('Uniform points')

plt.show()
