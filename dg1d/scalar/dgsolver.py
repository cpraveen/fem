import numpy as np
import matplotlib.pyplot as plt
from basis import *
from limiter import *

# SSPRK3 coefficients
ark = np.array([0.0, 3.0/4.0, 1.0/3.0])
brk = 1.0 - ark

class DG:
    def __init__(self,args):
        self.k = args.degree
        self.cfl = args.cfl / (2 * args.degree + 1)
        self.nc = args.ncell
        self.Tf = args.Tf
        self.plot_freq = args.plot_freq
        self.limit = args.limit
        self.compute_error = args.compute_error

        self.nd = args.degree + 1

        # Select initial condition
        if args.ic == 'sin2pi':
            from sin2pi import xmin,xmax,initial_condition
        elif args.ic == 'sin4pi':
            from sin4pi import xmin,xmax,initial_condition
        elif args.ic == 'gauss':
            from gauss import xmin,xmax,initial_condition
        elif args.ic == 'hat':
            from hat import xmin,xmax,initial_condition
        elif args.ic == 'mult':
            from mult import xmin,xmax,initial_condition
        else:
            print('Unknown initial condition')
            exit()

        self.xmin, self.xmax = xmin, xmax
        self.ic = initial_condition
        self.dx = (xmax - xmin) / self.nc
        self.Mdx2 = args.tvbM * self.dx**2

        # Select PDE
        if args.pde == 'linear':
            from linadv import flux,central_flux,upwind_flux,max_speed
            if args.num_flux == 'central':
                self.numflux = central_flux
            else:
                self.numflux = upwind_flux
        elif args.pde == 'varadv':
            from varadv import flux,numflux,max_speed
            self.numflux = numflux
        elif args.pde == 'burger':
            from burger import flux,central_flux,roe_flux,godunov_flux,max_speed
            if args.num_flux == 'central':
                self.numflux = central_flux
            elif args.num_flux == 'roe' or args.num_flux == 'upwind':
                self.numflux = roe_flux
            else:
                self.numflux = godunov_flux
        else:
            print('PDE not implemented')
            exit()

        self.flux = flux
        self.max_speed = max_speed

    def setup(self):
        nc, k, nd = self.nc, self.k, self.nd

        Nq = k + 1
        xg, wg = np.polynomial.legendre.leggauss(Nq)

        # Construct Vandermonde matrix for gauss points
        Vf = Vandermonde(k, xg)
        Vg = np.empty((Nq,nd))
        for i in range(Nq):
            for j in range(nd):
                Vg[i,j] = shape_grad (i, xg[j]) * wg[j]

        # Construct Vandermonde matrix for uniform points
        # uniform points in cell for plotting
        nu = np.max([2,k+1])
        xu = np.linspace(-1.0,+1.0,nu)
        Vu = Vandermonde(k, xu)

        # Required to evaluate solution at face
        bm, bp = np.empty(nd), np.empty(nd)
        for i in range(nd):
            bm[i] = shape_value(i,-1.0)
            bp[i] = shape_value(i,+1.0)

        self.Vf = Vf
        self.Vg = Vg
        self.Vu = Vu
        self.bm = bm
        self.bp = bp
        self.xg = xg
        self.wg = wg
        self.xu = xu

        # Allocate solution variables
        self.u0 = np.zeros((nc,nd)) # solution at n
        self.u1 = np.zeros((nc,nd)) # solution at n+1
        self.res= np.zeros((nc,nd)) # residual

    # Set initial condition by L2 projection
    def set_ic(self):
        nc = self.nc
        nd = self.nd
        dx = self.dx
        xg = self.xg
        wg = self.wg
        Vf = self.Vf
        ic = self.ic
        xmin = self.xmin

        for i in range(nc):
            xc = xmin + i*dx + 0.5*dx  # cell center
            x = xc + 0.5*dx*xg         # transform gauss points to real cell
            val = ic(x)
            for j in range(nd):
                self.u1[i, j] = 0.5 * val.dot(Vf[:, j]*wg)

    # Initialize plot
    def init_plot(self, ax):
        nc, xmin, xmax, dx = self.nc, self.xmin, self.xmax, self.dx
        xu, Vu = self.xu, self.Vu
        u1 = self.u1

        lines = []
        umin, umax = 1.0e20, -1.0e20
        for i in range(nc):  # Loop over cells
            xc = xmin + i*dx + 0.5*dx  # cell center
            x = xc + 0.5*dx*xu       # transform uniform points to cell
            f = Vu @ u1[i, :]
            line, = ax.plot(x, f, linewidth=2)
            lines.append(line)
            umin = np.min([umin, f.min()])
            umax = np.max([umax, f.max()])
        plt.title('Initial condition')
        plt.xlabel('x')
        plt.ylabel('u')
        plt.grid(True)
        plt.axis([xmin, xmax, umin-0.1, umax+0.1])
        plt.draw()
        plt.pause(0.1)
        self.lines = lines

    # Update plot
    def update_plot(self, t):
        nc, xmin, xmax, dx, Vu = self.nc, self.xmin, self.xmax, self.dx, self.Vu
        k, cfl = self.k, self.cfl
        lines = self.lines
        u1 = self.u1

        umin, umax = 1.0e20, -1.0e20
        for i in range(nc):
            xc = xmin + i*dx + 0.5*dx  # cell center
            f = Vu @ u1[i, :]
            lines[i].set_ydata(f)
            umin = np.min([umin, f.min()])
            umax = np.max([umax, f.max()])
        plt.axis([xmin, xmax, umin-0.1, umax+0.1])
        plt.title(str(nc)+' cells, Deg = '+str(k)+', CFL = '+str(round(cfl, 3)) +
                ', t = '+str("%.3f" % round(t, 3)))
        plt.draw()
        plt.pause(0.1)

    # Computes -rhs
    def compute_residual(self):
        nc = self.nc
        xmin, dx, xg, Vf, Vg = self.xmin, self.dx, self.xg, self.Vf, self.Vg
        bm, bp = self.bm, self.bp
        u1, res = self.u1, self.res
        flux = self.flux
        numflux = self.numflux

        # Loop over cells and compute cell integral
        for i in range(nc):
            xc = xmin + i*dx + 0.5*dx  # cell center
            x = xc + 0.5*dx*xg         # transform gauss points to cell
            u = Vf @ u1[i, :]          # solution at gauss points
            f = flux(x, u)             # flux at gauss points
            res[i, :] = -Vg @ f        # flux integral over cell

        # Now we compute the inter-cell fluxes
        # First face: left cell = last cell, right cell = 0'th cell
        ul = u1[-1, :].dot(bp)  # get ul from last cell
        ur = u1[ 0, :].dot(bm)  # get ur from 0'th cell
        f = numflux(xmin, ul, ur)
        res[-1, :] += f*bp  # Add to last cell
        res[ 0, :] -= f*bm  # Add to first cell

        # Loop over internal faces
        # Left cell = i-1, right cell = i
        for i in range(1, nc):
            ul = u1[i-1, :].dot(bp)
            ur = u1[i, :].dot(bm)
            f = numflux(xmin+i*dx, ul, ur)
            res[i-1, :] += f*bp
            res[i, :] -= f*bm

        # Last face is already done above using periodicity

    # Apply TVB limiter, using periodicity
    def apply_limiter(self):
        nc, limit, Mdx2 = self.nc, self.limit, self.Mdx2
        u1 = self.u1

        sqrt3 = np.sqrt(3.0)

        # Apply TVB limiter
        if limit == 'yes':
            for i in range(nc):
                if i == 0:
                    ul = u1[-1, 0]
                    ur = u1[1, 0]
                elif i == nc-1:
                    ul = u1[nc-2, 0]
                    ur = u1[0, 0]
                else:
                    ul = u1[i-1, 0]
                    ur = u1[i+1, 0]
                du = (1.0/sqrt3) * minmod(sqrt3*u1[i, 1], (u1[i, 0]-ul),
                                          (ur-u1[i, 0]), Mdx2)
                if np.abs(du-u1[i, 1]) > 1.0e-6:
                    u1[i, 1] = du   # Copy limited gradient
                    u1[i, 2:] = 0.0  # Kill all higher modes

    # Compute error norm using ng-point quadrature
    # We assume final solution = initial solution
    def compute_error_norm(self):
        k, nc, xmin, dx = self.k, self.nc, self.xmin, self.dx
        ic = self.ic
        u1 = self.u1

        ng = k + 3
        xg, wg = np.polynomial.legendre.leggauss(ng)
        Vf = Vandermonde(k, xg)

        error_norm = 0.0
        for i in range(nc):
            # DG solution at gauss points
            un = Vf @ u1[i, :]
            # Exact solution at gauss points
            xc = xmin + i*dx + 0.5*dx  # cell center
            x = xc + 0.5*dx*xg       # transform gauss points to cell
            ue = ic(x)
            error_norm += 0.5 * dx * np.sum((un-ue)**2 * wg)

        error_norm = np.sqrt(error_norm)
        return nc, dx, error_norm

    # Apply DG until final time
    def solve(self):
        cfl, dx, Tf = self.cfl, self.dx, self.Tf
        u0, u1, res = self.u0, self.u1, self.res
        max_speed = self.max_speed
        plot_freq = self.plot_freq
        limit = self.limit

        it, t = 0, 0.0
        while t < Tf:
            dt = cfl * dx / max_speed(u1)
            lam = dt / dx
            if t+dt > Tf:
                dt = Tf - t
                lam = dt / dx
            u0[:, :] = u1
            for rk in range(3):
                self.compute_residual()
                u1[:, :] = ark[rk]*u0 + brk[rk]*(u1 - lam*res)
                if limit == 'yes':
                    self.apply_limiter()
            t += dt; it += 1
            if (plot_freq > 0 and 
                (it % plot_freq == 0 or np.abs(Tf-t) < 1.0e-13)):
                self.update_plot(t)

    # Driver function
    def run(self):
        self.setup()
        self.set_ic()

        # plot initial condition
        if self.plot_freq > 0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            self.init_plot(ax)
            wait = input("Press enter to continue ")

        self.solve()
        if self.compute_error == 'yes':
            return self.compute_error_norm()
