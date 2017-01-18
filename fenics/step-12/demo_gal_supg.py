"""
Scalar advection on circles using Galerkin and SUPG
a.grad(u) = 0 in [0,1]x[0,1],  a = (y, -x)
u = uinlet on x=0
  = 0      on y=1
"""
from dolfin import *

def Inlet(x, on_boundary):
   return x[0] < DOLFIN_EPS and on_boundary

def Top(x, on_boundary):
   return x[1]-1 > -DOLFIN_EPS and on_boundary

class uinlet(Expression):
   def eval(self, value, x):
      if abs(x[1]-0.5) < 0.25:
         value[0]=1.0
      else:
         value[0]=0.0

np = 50

mesh = UnitSquare(np, np)
h = CellSize(mesh)

X = FunctionSpace(mesh, "CG", 1)

v = TestFunction(X)
u = TrialFunction(X)

# Note: moda = |a|
# If you change "a", then change "moda" also.
a    = Expression(("x[1]", "-x[0]"))
moda = Expression("sqrt(x[0]*x[0] + x[1]*x[1])")

a_gal = inner(a, grad(u))*v*dx
a_supg= (h/moda)*inner(a, grad(u))*inner(a, grad(v))*dx
#A     = a_gal + a_supg
A     = a_gal
L     = Constant(0)*v*dx

bc_in  = DirichletBC(X, uinlet(), Inlet)
bc_top = DirichletBC(X, 0, Top)
bc = [bc_in, bc_top]

u = Function(X)
solve(A == L, u, bc)
File("sol.pvd") << u

plot(u)
interactive()
