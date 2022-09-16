import numpy as np
from bucklev import *
import argparse
from scipy import optimize,integrate
from numpy import sqrt

parser = argparse.ArgumentParser()
parser.add_argument('-Tf', type = float, help = 'Time of plotting solution',
                    default = 0.4)
parser.add_argument('-ncell', type = int, help = 'Number of cells', default = 100)
args = parser.parse_args()

xmin, xmax = -1.0, 1.0
def hatbuck(x):
    f = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx < -0.5 or xx > 0.0:
            f[i] = 0.0
        else:
            f[i] = 1.0
    return f

# Finding the exact solution

u_s  = np.sqrt(a_buck/(1.0+a_buck)) # (f(u_s)-f(0))/(u_s-0)=f'(u_s)
u_ss = 1.0-1.0/np.sqrt(1.0+a_buck) # (f(u_ss)-f(1))/(u_ss-1)=f'(u_ss)

# Inverse of f' for computing rarefactions
# Inverse of f' restricted to [u_buck,1], an interval that contains [u_s,1]
def inv_f_s(v):
    # Inverse of f' at v equals root of this polynomial in [0.5,1]
    def p(u):
        value = v*(u**2+a_buck*(1.0-u)**2)**2-2.0*a_buck*u*(1.0-u)
        return value
    output = optimize.brentq(p, u_s, 1.0) # Gives root of polynomial
    return output

# Inverse of f' restricted to [0,u_buck], an interval that contains [0,u_ss]
def inv_f_ss(v):
    # Inverse of f' at v equals root of this polynomial in [0.5,1]
    def p(u):
        value = v*(u**2+a_buck*(1.0-u)**2)**2-2.0*a_buck*u*(1.0-u)
        return value
    output = optimize.brentq(p, 0.0, u_buck) # Gives root of polynomial
    return output

# To compute the exact solution, we need all shock locations for all time
# The shock location between the two rarefaction waves needs to be computed
# by solving the ode s'(t)=(f_l-f_r)/(u_l-u_r) for the particular time
def update_shock(shock,t,dt):
    def rh(shock,t):
      u_l = inv_f_ss((shock+0.5)/t)
      u_r = inv_f_s(shock/t)
      f_l = flux(0.0,u_l)
      f_r = flux(0.0,u_r)
      dsdt = (f_l-f_r)/(u_l-u_r)
      return dsdt
    # This is the time where rarefaction characteristics intersect
    if t>=1./(2.*fprime(u_ss)):
      time = [t,t+dt]
      output = integrate.odeint(rh,shock,time,rtol = 1e-5)
      return output[1]
    else:
      return 0.0

def exact_soln_hatbuck(x,t,shock):
  f = np.empty_like(x)
  f_u_s,f_u_ss = fprime(u_s),fprime(u_ss) # f'(u_s),f'(u_ss)
  for i,xx in enumerate(x):
    if xx <= -0.5:
      f[i] = 0.0
    elif xx>-0.5 and xx <= -0.5+f_u_ss*t:
      f[i] = inv_f_ss((xx-(-0.5))/t)
    elif xx>=-0.5+f_u_ss*t and xx <= 0.0:
      f[i] = 1.0
    elif xx>=0.0 and xx <=f_u_s*t:
      if xx>= shock:
        f[i] = inv_f_s(xx/t)
      else:
        f[i] = inv_f_ss((xx-(-0.5))/t)
    else:
      f[i] = 0.0
  return f

Tf = args.Tf
dt = Tf/1000.0
t = 0.0
shock = 0.0
while t < Tf:
    if (t+dt)>Tf:
        dt = Tf-t
    shock = update_shock(shock, t, dt) # Gives shock for solution at t+dt
    t += dt

M = args.ncell
dx = 2.0/M
xc = xmin + np.arange(M)*dx + 0.5*dx # cell centers
np.savetxt('exact.txt', np.column_stack([xc, exact_soln_hatbuck(xc, t, shock)]))
