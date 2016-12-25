# Plots P1 basis functions
import numpy as np
import matplotlib.pyplot as plt

def basis(x,i,h):
    xi = i*h; xl,xr = xi-h,xi+h
    y = np.empty_like(x)
    for j in range(len(x)):
        if x[j]>xr or x[j] < xl:
            y[j] = 0
        elif x[j]>xi:
            y[j] = (xr - x[j])/h
        else:
            y[j] = (x[j] - xl)/h
    return y

N = 5 # Number of elements
xg = np.linspace(0.0, 1.0, N+1)
h = 1.0/N

# Finer grid for plotting purpose
x = np.linspace(0.0,1.0,1000)

fig = plt.figure()
#plt.plot(xg,0*xg,'o-')
ax = fig.add_subplot(111)
ax.set_xlabel('$x$'); ax.set_ylabel('$\phi$')
line1, = ax.plot(xg,0*xg,'o-')
line2, = ax.plot(x,0*x,'r-',linewidth=2)
plt.axis([0.0,1.0,-0.1,1.5])
plt.draw()

for i in range(N+1):
    y = basis(x,i,h)
    line2.set_ydata(y)
    t = '$\phi_'+str(i)+'$'
    plt.title(t)
    plt.draw(); plt.pause(1.0)
