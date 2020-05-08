"""
-Laplace(u) = 1 in (0,1)x(0,1)
         u  = 0 on boundary
"""
from dolfin import *

mesh = UnitSquareMesh(20,20)

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Constant(1)
L = f*v*dx

# Dirichlet bc
g = Constant(0)
bc= DirichletBC(V, g, 'on_boundary')

# Solution variable
u = Function(V)
solve(a == L, u, bc)
u.rename('u','u')
File("sol.pvd") << u
