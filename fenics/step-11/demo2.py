"""
Heat equation in one space variable using Crank-Nicholson scheme
u_t = u_xx in (0, 2*pi) x (0, T)
u(x,0) = 0
u(0,t) = u(2*pi,t) = sin(t)
"""
from dolfin import *

n = 20
mesh = Interval(n, 0, 2*pi)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Solution variable
uold = Function(V)
unew = Function(V)

# Initial condition
uold.vector()[:] = 0.0

T = 4*pi # Final time
nt= 200  # Number of time steps
dt= T/nt # Time step

uavg = 0.5*(uold + u)
form = (1.0/dt)*(u - uold)*v*dx + inner(grad(uavg), grad(v))*dx
a    = lhs(form) # bilinear form
L    = rhs(form) # linear functional, zero in this case
A    = assemble(a)

# Time dependent bc at x=0 and x=2*pi
ubc  = Expression('sin(t)', t=0)
bc   = DirichletBC(V, ubc, "near(x[0],0) || near(x[0],2*pi)")

fsol = File('sol.pvd')
fsol << uold

t = 0
while t < T:
   b     = assemble(L)
   ubc.t = t + dt
   bc.apply(A, b)
   solve(A, unew.vector(), b)
   uold.assign(unew)
   t = t + dt
   print "t = ", t
   fsol << unew
