import numpy as np
import matplotlib.pyplot as plt
import argparse
from basis import *
from limiter import *

# SSPRK3 coefficients
ark = np.array([0.0, 3.0/4.0, 1.0/3.0])
brk = 1.0 - ark

sqrt3  = np.sqrt(3.0)
isqrt3 = 1.0/sqrt3

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ncell', type=int, help='Number of cells', default=100)
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-Tf', type=float, help='Final time', default=0.0)
parser.add_argument('-plot_freq', type=int, help='Frequency to plot solution', 
                    default=5)
parser.add_argument('-ic', choices=('sod','lax'), help='Initial condition', 
                    default='sod')
parser.add_argument('-char_lim', type=int, help='Characteristic limiter', 
                    default=0)
parser.add_argument('-tvbM', type=float, help='TVB M parameter', default=0.0)
args = parser.parse_args()

# Select initial condition
if args.ic == 'sod':
    from sod import *
elif args.ic == 'lax':
    from lax import *
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
def init_plot(fig,ax,rho,mom,ene):
    lines0,lines1,lines2 = [], [], []
    umin0, umax0 = 1.0e20, -1.0e20
    umin1, umax1 = 1.0e20, -1.0e20
    umin2, umax2 = 1.0e20, -1.0e20
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        # Density
        r  = Vu.dot(rho[i,:])
        line, = ax[0].plot(x,r,linewidth=2)
        lines0.append(line)
        umin0 = np.min([umin0, r.min()])
        umax0 = np.max([umax0, r.max()])
        # Velocity
        m  = Vu.dot(mom[i,:])
        v  = m / r
        line, = ax[1].plot(x,v,linewidth=2)
        lines1.append(line)
        umin1 = np.min([umin1, v.min()])
        umax1 = np.max([umax1, v.max()])
        # Pressure
        E  = Vu.dot(ene[i,:])
        pre= pressure(r, m, E)
        line, = ax[2].plot(x,pre,linewidth=2)
        lines2.append(line)
        umin2 = np.min([umin2, pre.min()])
        umax2 = np.max([umax2, pre.max()])
    fig.canvas.manager.set_window_title('Initial condition')
    ax[0].set_xlabel('x'); ax[0].set_title('Density'); ax[0].grid(True)
    ax[0].set_xlim([xmin,xmax])
    ax[0].set_ylim([umin0-0.05,umax0+0.05])
    ax[1].set_xlabel('x'); ax[1].set_title('Velocity'); ax[1].grid(True)
    ax[1].set_xlim([xmin,xmax])
    ax[1].set_ylim([umin1-0.05,umax1+0.05])
    ax[2].set_xlabel('x'); ax[2].set_title('Pressure'); ax[2].grid(True)
    ax[2].set_xlim([xmin,xmax])
    ax[2].set_ylim([umin2-0.05,umax2+0.05])
    plt.draw(); plt.pause(0.1)
    return lines0,lines1,lines2

# Update plot
def update_plot(fig,ax,lines0,lines1,lines2,t,rho,mom,ene):
    umin0, umax0 = 1.0e20, -1.0e20
    umin1, umax1 = 1.0e20, -1.0e20
    umin2, umax2 = 1.0e20, -1.0e20
    for i in range(nc):
        xc = xmin + i*dx + 0.5*dx # cell center
        x  = xc + 0.5*dx*xu       # transform gauss points to cell
        # Density
        r  = Vu.dot(rho[i,:])
        lines0[i].set_ydata(r)
        umin0 = np.min([umin0, r.min()])
        umax0 = np.max([umax0, r.max()])
        # Velocity
        m  = Vu.dot(mom[i,:])
        v  = m / r
        lines1[i].set_ydata(v)
        umin1 = np.min([umin1, v.min()])
        umax1 = np.max([umax1, v.max()])
        # Pressure
        E  = Vu.dot(ene[i,:])
        pre= pressure(r, m, E)
        lines2[i].set_ydata(pre)
        umin2 = np.min([umin2, pre.min()])
        umax2 = np.max([umax2, pre.max()])
    fig.canvas.manager.set_window_title(str(nc)+' cells, Deg = '+str(k)+
                                        ', CFL = '+str(cfl)+', t = '+('%.3e'%t))
    ax[0].set_ylim([umin0-0.05,umax0+0.05])
    ax[1].set_ylim([umin1-0.05,umax1+0.05])
    ax[2].set_ylim([umin2-0.05,umax2+0.05])
    plt.draw(); plt.pause(0.1)

# Allocate solution variables
rho0 = np.zeros((nc,nd)) # solution at n
rho1 = np.zeros((nc,nd)) # solution at n+1
mom0 = np.zeros((nc,nd)) # solution at n
mom1 = np.zeros((nc,nd)) # solution at n+1
ene0 = np.zeros((nc,nd)) # solution at n
ene1 = np.zeros((nc,nd)) # solution at n+1
resr = np.zeros((nc,nd)) # mass residual
resm = np.zeros((nc,nd)) # momentum residual
rese = np.zeros((nc,nd)) # energy residual

# Set initial condition by L2 projection
for i in range(nc):
    xc = xmin + i*dx + 0.5*dx # cell center
    x  = xc + 0.5*dx*xg       # transform gauss points to cell
    rho,mom,ene = initial_condition(x)
    for j in range(nd):
        rho1[i,j] = 0.5 * rho.dot(Vf[:,j]*wg)
        mom1[i,j] = 0.5 * mom.dot(Vf[:,j]*wg)
        ene1[i,j] = 0.5 * ene.dot(Vf[:,j]*wg)

# plot initial condition
fig,ax = plt.subplots(nrows=1,ncols=3,figsize=(15,5))
lines0,lines1,lines2 = init_plot(fig,ax,rho1,mom1,ene1)
wait = input("Press enter to continue ")

# If final time given on command line, use that
if args.Tf > 0.0:
    Tf  = args.Tf

it, t = 0, 0.0
while t < Tf:
    dt  = cfl*dx/max_speed(rho1,mom1,ene1)
    lam = dt/dx
    if t+dt > Tf:
        dt = Tf - t
        lam = dt/dx
    rho0[:,:] = rho1
    mom0[:,:] = mom1
    ene0[:,:] = ene1
    for rk in range(3):
        # Loop over cells and compute cell integral
        for i in range(nc):
            xc = xmin + i*dx + 0.5*dx # cell center
            x  = xc + 0.5*dx*xg       # transform gauss points to cell
            rho = Vf.dot(rho1[i,:])   # solution at gauss points
            mom = Vf.dot(mom1[i,:])   # solution at gauss points
            ene = Vf.dot(ene1[i,:])   # solution at gauss points
            frho,fmom,fene = flux(rho,mom,ene) # flux at gauss points
            resr[i,:] = -Vg.dot(frho)
            resm[i,:] = -Vg.dot(fmom)
            rese[i,:] = -Vg.dot(fene)
        # First face: neumann bc
        rhol = rhor = rho1[0,:].dot(bm)
        moml = momr = mom1[0,:].dot(bm)
        enel = ener = ene1[0,:].dot(bm)
        f  = numflux([rhol,moml,enel],[rhor,momr,ener])
        resr[0,:] -= f[0]*bm # Add to first cell
        resm[0,:] -= f[1]*bm # Add to first cell
        rese[0,:] -= f[2]*bm # Add to first cell
        # Loop over internal faces
        # Left cell = i-1, right cell = i
        for i in range(1,nc):
            rhol = rho1[i-1,:].dot(bp)
            rhor = rho1[i  ,:].dot(bm)
            moml = mom1[i-1,:].dot(bp)
            momr = mom1[i  ,:].dot(bm)
            enel = ene1[i-1,:].dot(bp)
            ener = ene1[i  ,:].dot(bm)
            f  = numflux([rhol,moml,enel], [rhor,momr,ener])
            resr[i-1,:] += f[0]*bp # to left cell
            resr[i  ,:] -= f[0]*bm # to right cell
            resm[i-1,:] += f[1]*bp # to left cell
            resm[i  ,:] -= f[1]*bm # to right cell
            rese[i-1,:] += f[2]*bp # to left cell
            rese[i  ,:] -= f[2]*bm # to right cell
        # last face: neumann bc
        rhol = rhor = rho1[-1,:].dot(bp)
        moml = momr = mom1[-1,:].dot(bp)
        enel = ener = ene1[-1,:].dot(bp)
        f  = numflux([rhol,moml,enel],[rhor,momr,ener])
        resr[-1,:] += f[0]*bp # Add to last cell
        resm[-1,:] += f[1]*bp # Add to last cell
        rese[-1,:] += f[2]*bp # Add to last cell
        # Peform rk stage
        rho1[:,:] = ark[rk]*rho0 + brk[rk]*(rho1 - lam*resr)
        mom1[:,:] = ark[rk]*mom0 + brk[rk]*(mom1 - lam*resm)
        ene1[:,:] = ark[rk]*ene0 + brk[rk]*(ene1 - lam*rese)
        # Apply limiter
        if k > 0:
            dul,uc,dur,du,dun = np.zeros(3),np.zeros(3),np.zeros(3),\
                                np.zeros(3),np.zeros(3)
            for i in range(1,nc-1):
                dul[0] = rho1[i  ,0] - rho1[i-1,0] # density
                dur[0] = rho1[i+1,0] - rho1[i  ,0]
                uc[0]  = rho1[i  ,0]
                du[0]  = rho1[i  ,1]
                dul[1] = mom1[i  ,0] - mom1[i-1,0] # momentum
                dur[1] = mom1[i+1,0] - mom1[i  ,0]
                uc[1]  = mom1[i  ,0]
                du[1]  = mom1[i  ,1]
                dul[2] = ene1[i  ,0] - ene1[i-1,0] # energy
                dur[2] = ene1[i+1,0] - ene1[i  ,0]
                uc[2]  = ene1[i  ,0]
                du[2]  = ene1[i  ,1]
                # Convert to characteristic variables
                if args.char_lim == 1:
                    R, L = EigMatrix(uc)
                    du   = L.dot(du)
                    dul  = L.dot(dul)
                    dur  = L.dot(dur)
                # Check if any slope has changed
                lim   = 0
                for j in range(3):
                    dun[j] = isqrt3*minmod(sqrt3*du[j],dul[j],dur[j],Mdx2)
                    if np.abs(dun[j]-du[j]) > 1.0e-6:
                        lim = 1
                if lim == 1:
                    # Convert back to conserved
                    if args.char_lim == 1:
                        dun = R.dot(dun)
                    rho1[i,1] = dun[0]; rho1[i,2:] = 0.0
                    mom1[i,1] = dun[1]; mom1[i,2:] = 0.0
                    ene1[i,1] = dun[2]; ene1[i,2:] = 0.0
    t += dt; it += 1
    if it%args.plot_freq == 0 or np.abs(Tf-t) < 1.0e-13:
        update_plot(fig,ax,lines0,lines1,lines2,t,rho1,mom1,ene1)

plt.show() # Dont close window at end of program
