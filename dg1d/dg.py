import numpy as np
import matplotlib.pyplot as plt
from basis import *
from linadv import *

# SSPRK3 coefficients
ark = np.array([0.0, 3.0/4.0, 1.0/3.0])
brk = 1.0 - ark

def initial_condition(x):
    return np.sin(2*np.pi*x)

xmin, xmax = 0.0, 1.0

cfl= 0.8
k  = 1     # polynomial degree
nd = k + 1 # dofs per cell
nc = 50    # number of cells
dx = (xmax - xmin)/nc

# k+1 gauss rule, integrate exactly upto degree 2*k+1
xg, wg = np.polynomial.legendre.leggauss(k+1)
wg *= 0.5

# Construct Vandermonde matrix for gauss points
Vf = np.zeros((nd,nd))
Vg = np.zeros((nd,nd))
for i in range(nd):
    for j in range(nd):
        Vf[i,j] = shape_value(j, xg[i])
        Vg[i,j] = shape_grad (j, xg[i])

# Construct Vandermonde matrix for uniform points
# uniform points in cell for plotting
nu = np.max([2,k+1])
xu = np.linspace(-1.0,+1.0,nu)
Vu = np.zeros((nu,k+1))
for i in range(nu):
    for j in range(k+1):
        Vu[i,j] = shape_value(j, xu[i])

u0 = np.zeros((nc,nd)) # solution at n
u1 = np.zeros((nc,nd)) # solution at n+1
res= np.zeros((nc,nd)) # solution at n+1

# Set initial condition by L2 projection
for i in range(nc):
    xc = xmin + i*dx + 0.5*dx # cell center
    x  = xc + 0.5*dx*xg       # transform gauss points to cell
    val= initial_condition(x)
    for j in range(nd):
        u0[i,j] = val.dot(Vf[:,j]*wg)

u1[:,:] = u0

def init_plot(ax,u0):
    lines = []
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu.dot(u0[i,:])
        line, = ax.plot(x,f)
        lines.append(line)
    return lines

def update_plot(lines,t,u1):
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu.dot(u1[i,:])
        lines[i].set_ydata(f)
    plt.title('t = '+('%.3e'%t))

# plot initial condition
fig = plt.figure()
ax = fig.add_subplot(111)
lines = init_plot(ax,u0)
plt.title('Initial condition')
plt.xlabel('x')
plt.ylabel('u')
plt.grid(True)
plt.draw(); plt.pause(0.1)
wait = raw_input("Press enter to continue ")

t = 0.0
dt= cfl*dx/(2*k+1)
Tf= 1.0
lam = dt/dx
while t < Tf:
    for rk in range(3):
        # Loop over cells and compute cell integral
        for i in range(nc):
            u = Vf.dot(u1[i,:]) # solution at gauss points
            f = flux(u)        # flux at gauss points
            for j in range(nd):
                res[i,j] = -2.0*f.dot(Vg[:,j]*wg)
        # First face: left cell = nc-1, right cell = 0
        ul = u1[-1,:].dot(Vu[-1,:])
        ur = u1[ 0,:].dot(Vu[ 0,:])
        f  = numflux(ul, ur)
        res[-1,:] += f*Vu[-1,:]
        res[ 0,:] -= f*Vu[ 0,:]
        # Loop over internal faces
        # Left cell = i-1, right cell = i
        for i in range(1,nc):
            ul = u1[i-1,:].dot(Vu[-1,:])
            ur = u1[i  ,:].dot(Vu[ 0,:])
            f  = numflux(ul, ur)
            res[i-1,:] += f*Vu[-1,:]
            res[i  ,:] -= f*Vu[ 0,:]
        # Peform rk stage
        u1[:,:] = ark[rk]*u0 + brk[rk]*(u1 - lam*res)
    u0[:,:] = u1
    t += dt
    update_plot(lines,t,u1); plt.draw(); plt.pause(0.1)
