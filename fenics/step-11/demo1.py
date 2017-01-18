"""
Heat equation in one space variable using Crank-Nicholson scheme
u_t = u_xx + f in (0, 2*pi) x (0, T)
f = sin(x)*(cos(t) - sin(t))
u(x,0) = sin(x)
u(0,t) = u(2*pi,t) = 0
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
uinit= Expression('sin(x[0])')
uold = interpolate(uinit, V)

T = 4*pi # Final time
nt=100   # Number of time steps
dt= T/nt # Time step

f = Expression('sin(x[0])*(cos(t) - sin(t))', t=0)
uavg = 0.5*(uold + u)
form = (1.0/dt)*(u - uold)*v*dx + inner(grad(uavg), grad(v))*dx - f*v*dx
a    = lhs(form)
L    = rhs(form)
A    = assemble(a)

# Homogeneous bc at x=0 and x=2*pi
bc   = DirichletBC(V, 0.0, "x[0] < DOLFIN_EPS || x[0]-2*pi > -DOLFIN_EPS")

fsol = File('sol.pvd')
fsol << uold

t = 0
while t < T:
   f.t = t + 0.5*dt
   b = assemble(L)
   bc.apply(A, b)
   solve(A, unew.vector(), b)
   uold.assign(unew)
   t = t + dt
   print "t = ", t
   fsol << unew
