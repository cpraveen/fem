"""
-Laplace(u) = 1 in (0,1)x(0,1)
         u  = y(1-y) on x=0
         u  = 0      on y=0 and y=1
      du/dn = 1      on x=1
"""
from dolfin import *

# Characteristic function for different boundaries
class Left(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0.0)

class Right(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 1.0)

class Bottom(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0)

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 1.0)

mesh = UnitSquareMesh(20,20)

boundaries = FacetFunction("uint", mesh)
boundaries.set_all(0)
right = Right()
right.mark(boundaries, 1)
dsn = Measure("ds")[boundaries]

V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Constant(1.0)
g = Constant(1.0)
L = f*v*dx + g*v*dsn(1)

# Dirichlet bc
u_left   = Expression('x[1]*(1.0-x[1])')
u_bottom = Constant(0.0)
u_top    = Constant(0.0)

bc_left   = DirichletBC(V, u_left,   Left())
bc_bottom = DirichletBC(V, u_bottom, Bottom())
bc_top    = DirichletBC(V, u_top,    Top())

bc = [bc_left, bc_bottom, bc_top]

# Solution variable
w = Function(V)

solve(a == L, w, bc)

file = File('sol.pvd')
file << w

plot(mesh)
plot(w)
interactive()
