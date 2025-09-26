"""
Solve u_t = u_xx in [0,1]
u(0,t) = u(1,t) = 0

See Fig 7 in:
Alhawary and Zhang, A study of DG methods for diffusion using the combined-mode analysis
"""
import numpy as np
import matplotlib.pyplot as plt
from firedrake import *

k    = 2        # degree of polynomial
Cip  = (k+1)**2 # IP penalty
nc   = 6        # number of points on each side of square

mesh = UnitIntervalMesh(nc)
V    = FunctionSpace(mesh, 'DG', k)
u    = TrialFunction(V)
v    = TestFunction(V)
n    = FacetNormal(mesh)
h    = CellDiameter(mesh)
havg = (h('-')+h('+'))/2.0

f = Constant(0.0) # RHS function
g = Constant(0.0) # Dirichlet bc

F =   inner(grad(u),grad(v))*dx         \
    - inner(avg(grad(u)), jump(v,n))*dS \
    - inner(avg(grad(v)), jump(u,n))*dS \
    - inner(grad(u), v*n)*ds            \
    - inner(grad(v), (u-g)*n)*ds        \
    + (Cip/havg)*jump(u)*jump(v)*dS \
    + (Cip/h)*(u-g)*v*ds            \
    - f*v*dx

uinit = Expression('sin(6*pi*x[0])',degree=2*k+1)
u0 = project(uinit, V)

Vf = FunctionSpace(mesh, 'DG', 10) # for visualization
uf = Functio(Vf).interpolate(u0)
ua=uf.vector().get_local()
x = SpatialCoordinate(mesh)
xf= Function(Vf).interpolate(x)
idx = np.argsort(xf[:,0])
plt.plot(xf[idx,0],ua[idx]);
plt.savefig('sol0.pdf')
print('See file sol0.pdf')

dt = 1.0e-4
Tf = (2.0/9.0)/nc**2
t = 0.0
nt = int(Tf/dt)
idt = Constant(1/dt)

# First step: BDF1
bdf1 = idt*(u - u0)*v*dx + F
a,L = lhs(bdf1), rhs(bdf1)
u1 = Function(V)
solve(a == L, u1)
t = t + dt

# From now on: BDF2
bdf2 = idt*(1.5*u - 2.0*u1 + 0.5*u0)*v*dx + F
a,L = lhs(bdf2), rhs(bdf2)
u2 = Function(V)

for i in range(nt):
    solve(a==L, u2)
    t = t + dt
    print('i, t =', i, t)
    u0.assign(u1)
    u1.assign(u2)

uf = Function(Vf).interpolate(u2)
ua=uf.vector().get_local()
plt.figure()
plt.plot(xf[idx,0],ua[idx]);
plt.savefig('sol.pdf')
print('See file sol.pdf')
