import numpy as np
import matplotlib.pyplot as plt
import argparse
from basis import *

# SSPRK3 coefficients
ark = np.array([0.0, 3.0/4.0, 1.0/3.0])
brk = 1.0 - ark

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pde', choices=('linear','burger'), help='PDE', 
                    default='linear')
parser.add_argument('-ncell', type=int, help='Number of cells', default=50)
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-save_freq', type=int, help='Frequency to save solution', 
                    default=1)
parser.add_argument('-ic', choices=('sine','hat'), help='Initial condition', 
                    default='sine')
parser.add_argument('-limit', choices=('no','yes'), help='Apply limiter', 
                    default='no')
args = parser.parse_args()

if args.pde == 'linear':
    from linadv import *
elif args.pde == 'burger':
    from burger import *
else:
    print('PDE not implemented')
    exit()

if args.ic == 'sine':
    from sine import *
elif args.ic == 'hat':
    from hat import *
else:
    print('Unknown initial condition')
    exit()

cfl= args.cfl     # cfl number
k  = args.degree  # polynomial degree
nc = args.ncell   # number of cells

nd = k + 1 # dofs per cell
dx = (xmax - xmin)/nc

# k+1 point gauss rule, integrates exactly upto degree 2*k+1
xg, wg = np.polynomial.legendre.leggauss(k+1)
wg *= 0.5 # normal weights so their sum is one

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
    umin, umax = 1.0e20, -1.0e20
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu.dot(u0[i,:])
        line, = ax.plot(x,f)
        lines.append(line)
        umin = np.min([umin, f.min()])
        umax = np.max([umax, f.max()])
    plt.title('Initial condition')
    plt.xlabel('x'); plt.ylabel('u'); plt.grid(True)
    plt.axis([xmin,xmax,umin-0.1,umax+0.1])
    plt.draw(); plt.pause(0.1)
    return lines

def update_plot(lines,t,u1):
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu.dot(u1[i,:])
        lines[i].set_ydata(f)
    plt.title(str(nc)+' cells, CFL = '+str(cfl)+', t = '+('%.3e'%t))
    plt.draw(); plt.pause(0.1)

# plot initial condition
fig = plt.figure()
ax = fig.add_subplot(111)
lines = init_plot(ax,u0)
wait = raw_input("Press enter to continue ")

it, t = 0, 0.0
dt  = cfl*dx/(2*k+1)/max_speed(u0)
Tf  = args.Tf
lam = dt/dx
while t < Tf:
    if t+dt > Tf:
        dt = Tf - f
        lam = dt/dx
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
        res[-1,:] += f*Vu[-1,:] # Add to last cell
        res[ 0,:] -= f*Vu[ 0,:] # Add to first cell
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
        if args.limit == 'yes':
            print('Limiter not implemented')
            exit()
    u0[:,:] = u1
    t += dt; it += 1
    if it%args.save_freq == 0:
        update_plot(lines,t,u1) 
