import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

# See section on "Gauss-Lobattto rules" here
#    https://en.wikipedia.org/wiki/Gaussian_quadrature
def gauss_lobatto_points(n):
    if n == 2:
        x = [-1,1]
    elif n == 3:
        x = [-1,0,1]
    elif n == 4:
        x = [-1, -1/np.sqrt(5), 1/np.sqrt(5), 1]
    elif n == 5:
        x = [-1, -np.sqrt(3.0/7.0), 0, np.sqrt(3.0/7.0), 1]
    elif n == 6:
        x = [-1, -np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
                 -np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
              1]
    elif n == 7:
        x = [-1, -np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
                 -np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  0,
                  np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
              1]
    else:
        print('Not implemented')
        exit()
    return np.array(x)

# Plot the points
for n in range(2,8):
    x = 0.5*(1 + gauss_lobatto_points(n))
    plt.plot(x, (10-n)*np.ones(n), 'o-')
plt.yticks([])
plt.title('Gauss-Lobatto points')

for n in range(2,8):
    plt.figure()
    xs = 0.5*(1 + gauss_lobatto_points(n))
    for i in range(n):
        values = np.zeros((n,1))
        values[i] = 1.0
        fun = lagrange(xs, values)
        x = np.linspace(0,1,100)
        plt.plot(x, fun(x))
        plt.xlabel('x')
        plt.grid(True)
        plt.title('Number of points = '+str(n))

plt.show()
