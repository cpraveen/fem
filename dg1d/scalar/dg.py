"""
Solve scalar conservation law with periodic bc
To get help, type
    python dg.py -h
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from basis import *
from limiter import *

# SSPRK3 coefficients
ark = np.array([0.0, 3.0/4.0, 1.0/3.0])
brk = 1.0 - ark

sqrt3 = np.sqrt(3.0)

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pde', choices=('linear','varadv','burger'), help='PDE', 
                    default='linear')
parser.add_argument('-ncell', type=int, help='Number of cells', default=50)
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution', 
                    default=1)
parser.add_argument('-ic', choices=('sin2pi','sin4pi','gauss','hat','mult'),
                    help='Initial condition', default='sin2pi')
parser.add_argument('-limit', choices=('no','yes'), help='Apply limiter', 
                    default='no')
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
parser.add_argument('-compute_error', choices=('no','yes'), 
                    help='Compute error norm', default='no')
parser.add_argument('-num_flux', choices=('central','upwind','roe','godunov'),
                    help='Numerical flux', default='upwind')
args = parser.parse_args()

# Select PDE
if args.pde == 'linear':
    from linadv import *
    if args.num_flux == 'central':
        numflux = central_flux
    else:
        numflux = upwind_flux
elif args.pde == 'varadv':
    from varadv import *
elif args.pde == 'burger':
    from burger import *
    if args.num_flux == 'central':
        numflux = central_flux
    elif args.num_flux == 'roe' or args.num_flux == 'upwind':
        numflux = roe_flux
    else:
        numflux = godunov_flux
else:
    print('PDE not implemented')
    exit()

# Select initial condition
if args.ic == 'sin2pi':
    from sin2pi import *
elif args.ic == 'sin4pi':
    from sin4pi import *
elif args.ic == 'gauss':
    from gauss import *
elif args.ic == 'hat':
    from hat import *
elif args.ic == 'mult':
    from mult import *
else:
    print('Unknown initial condition')
    exit()

k  = args.degree      # polynomial degree
cfl= args.cfl/(2*k+1) # cfl number
nc = args.ncell       # number of cells

nd = k + 1 # dofs per cell
dx = (xmax - xmin)/nc
Mdx2 = args.tvbM * dx**2

# k+1 point gauss rule, integrates exactly upto degree 2*k+1
Nq     = k+1
xg, wg = np.polynomial.legendre.leggauss(Nq)

# Construct Vandermonde matrix for gauss points
Vf = np.zeros((Nq,nd))
Vg = np.zeros((Nq,nd))
for i in range(Nq):
    for j in range(nd):
        Vf[i,j] = shape_value(j, xg[i])
        Vg[i,j] = shape_grad (i, xg[j]) * wg[j]

# Construct Vandermonde matrix for uniform points
# uniform points in cell for plotting
nu = np.max([2,k+1])
xu = np.linspace(-1.0,+1.0,nu)
Vu = np.zeros((nu,k+1))
for i in range(nu):
    for j in range(k+1):
        Vu[i,j] = shape_value(j, xu[i])

# Required to evaluate solution at face
bm, bp = np.zeros(nd), np.zeros(nd)
for i in range(nd):
    bm[i] = shape_value(i,-1.0)
    bp[i] = shape_value(i,+1.0)

# Initialize plot
def init_plot(ax,u0):
    lines = []
    umin, umax = 1.0e20, -1.0e20
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu @ u0[i,:]
        line, = ax.plot(x,f,linewidth=2)
        lines.append(line)
        umin = np.min([umin, f.min()])
        umax = np.max([umax, f.max()])
    plt.title('Initial condition')
    plt.xlabel('x'); plt.ylabel('u'); plt.grid(True)
    plt.axis([xmin,xmax,umin-0.1,umax+0.1])
    plt.draw(); plt.pause(0.1)
    return lines

# Update plot
def update_plot(lines,t,u1):
    umin, umax = 1.0e20, -1.0e20
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        f  = Vu @ u1[i,:]
        lines[i].set_ydata(f)
        umin = np.min([umin, f.min()])
        umax = np.max([umax, f.max()])
    plt.axis([xmin,xmax,umin-0.1,umax+0.1])
    plt.title(str(nc)+' cells, Deg = '+str(k)+', CFL = '+str(round(cfl,3))+
              ', t = '+str(round(t,3)))
    plt.draw(); plt.pause(0.1)

# Allocate solution variables
u0 = np.zeros((nc,nd)) # solution at n
u1 = np.zeros((nc,nd)) # solution at n+1
res= np.zeros((nc,nd)) # residual

# Set initial condition by L2 projection
for i in range(nc):
    xc = xmin + i*dx + 0.5*dx # cell center
    x  = xc + 0.5*dx*xg       # transform gauss points to real cell
    val= initial_condition(x)
    for j in range(nd):
        u1[i,j] = 0.5 * val.dot(Vf[:,j]*wg)

# plot initial condition
fig = plt.figure()
ax = fig.add_subplot(111)
lines = init_plot(ax,u1)
wait = input("Press enter to continue ")

it, t = 0, 0.0
dt  = cfl*dx/max_speed(u1)
Tf  = args.Tf
lam = dt/dx
while t < Tf:
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/dx
    u0[:,:] = u1
    for rk in range(3):
        # Loop over cells and compute cell integral
        for i in range(nc):
            xc = xmin + i*dx + 0.5*dx # cell center
            x  = xc + 0.5*dx*xg       # transform gauss points to cell
            u = Vf @ u1[i,:]          # solution at gauss points
            f = flux(x,u)             # flux at gauss points
            res[i,:] = -Vg @ f        # flux integral over cell
        # Now we compute the inter-cell fluxes
        # First face: left cell = last cell, right cell = 0'th cell
        ul = u1[-1,:].dot(bp) # get ul from last cell
        ur = u1[ 0,:].dot(bm) # get ur from 0'th cell
        f  = numflux(xmin, ul, ur)
        res[-1,:] += f*bp # Add to last cell
        res[ 0,:] -= f*bm # Add to first cell
        # Loop over internal faces
        # Left cell = i-1, right cell = i
        for i in range(1,nc):
            ul = u1[i-1,:].dot(bp)
            ur = u1[i  ,:].dot(bm)
            f  = numflux(xmin+i*dx, ul, ur)
            res[i-1,:] += f*bp
            res[i  ,:] -= f*bm
        # Peform rk stage
        u1[:,:] = ark[rk]*u0 + brk[rk]*(u1 - lam*res)
        # Apply TVB limiter
        if args.limit == 'yes':
            for i in range(nc):
                if i==0:
                    ul = u1[-1,0]
                    ur = u1[ 1,0]
                elif i==nc-1:
                    ul = u1[nc-2,0]
                    ur = u1[0   ,0]
                else:
                    ul = u1[i-1,0]
                    ur = u1[i+1,0]
                du = (1.0/sqrt3)*minmod(sqrt3*u1[i,1], (u1[i,0]-ul),
                                        (ur-u1[i,0]), Mdx2)
                if np.abs(du-u1[i,1]) > 1.0e-6:
                    u1[i,1 ] = du   # Copy limited gradient
                    u1[i,2:] = 0.0  # Kill all higher modes
    t += dt; it += 1
    if it%args.plot_freq == 0 or np.abs(Tf-t) < 1.0e-13:
        update_plot(lines,t,u1)

if args.compute_error == 'yes':
    # Compute error norm using ng-point quadrature
    # We assume final solution = initial solution
    ng = k + 3
    xg, wg = np.polynomial.legendre.leggauss(ng)

    Vf = np.zeros((ng,k+1))
    for i in range(ng):
        for j in range(k+1):
            Vf[i,j] = shape_value(j, xg[i])

    error_norm = 0.0
    for i in range(nc):
        # DG solution at gauss points
        un = Vf @ u1[i,:]
        # Exact solution at gauss points
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xg       # transform gauss points to cell
        ue = initial_condition(x)
        error_norm += 0.5 * dx * np.sum( (un-ue)**2 * wg )

    print('L2 error norm = %e\n' % np.sqrt(error_norm))

plt.show() # Dont close window at end of program
