# Plots P2 basis functions
import numpy as np
import matplotlib.pyplot as plt

def basis1(h,x):
    y = np.empty_like(x)
    for i in range(len(x)):
        if x[i] < 0.0:
            y[i] = 2*(x[i]+h)*(x[i]+h/2)/h**2
        else:
            y[i] = 2.0*(h/2-x[i])*(h-x[i])/h**2
    return y

def basis2(h,x):
    return 4*x*(h-x)/h**2

N = 5 # Number of elements
xg = np.linspace(0.0, 1.0, N+1)
h = 1.0/N

M = 2*N + 1 # no of dofs
xd = np.linspace(0.0, 1.0, M)


x1 = np.linspace(-h,h,500)
x2 = np.linspace(0.0,h,500)
y1 = basis1(h,x1)
y2 = basis2(h,x2)

xs1,xs2 = 0.0,0.0
for i in range(M):
    plt.clf()
    plt.plot(xg,0*xg,'s-',markersize=10)
    plt.plot(xd,0*xd,'ro',markersize=8)
    plt.axis([0.0,1.0,-0.5,1.5])
    plt.xlabel('$x$'); plt.ylabel('$\phi$')
    if i%2==0:
        plt.plot(x1+xs1,y1,linewidth=2)
        xs1 += h
    else:
        plt.plot(x2+xs2,y2,linewidth=2)
        xs2 += h
    # Print global dof numbers
    for j in range(M):
        plt.text(xd[j],-0.15,str(j))
    t = '$\phi_{'+str(i)+'}$'
    plt.title(t)
    plt.draw(); plt.pause(2.0)
