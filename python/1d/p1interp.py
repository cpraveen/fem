import numpy as np
import matplotlib.pyplot as plt

n = 20
xmin, xmax = 0.0, 1.0

f = lambda x: 2.0 + np.sin(15*x) * np.exp(-x)

x = np.linspace(xmin, xmax, 10*n)
xg = np.linspace(xmin, xmax, n)

plt.figure()
plt.plot(x, f(x), '-')
plt.ylim(-0.1,3)
plt.xlabel('x'); plt.ylabel('f(x)'); plt.grid(True)
plt.title('Some function $f(x)$')
plt.show()

plt.figure()
plt.plot(x, f(x), '-')
plt.plot(xg, 0*xg, 'o')
plt.ylim(-0.10,3)
plt.xlabel('x'); plt.ylabel('f(x)'); plt.grid(True)
plt.title('Partition domain with disjoint elements')
plt.show()

plt.figure()
for i in range(n):
    plt.plot([xg[i],xg[i]],[0,f(xg[i])],'k-')
plt.plot(x, f(x), '-')
plt.plot(xg, f(xg), 's')
plt.plot(xg, 0*xg, 'o')
plt.ylim(-0.10,3)
plt.xlabel('x'); plt.ylabel('f(x)'); plt.grid(True)
plt.title('Evaluate function on the grid')
plt.show()

plt.figure()
plt.plot(xg, f(xg), '-s')
plt.plot(xg, 0*xg, 'o')
plt.ylim(-0.10,3)
plt.xlabel('x'); plt.ylabel('f(x)'); plt.grid(True)
plt.title('Approx $f(x)$ by piecewise linear function')
plt.show()
