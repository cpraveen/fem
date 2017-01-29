"""
-Laplace(u) = 1 in (0,1)x(0,1)
         u  = 0 on boundary
"""
from dolfin import *

def Boundary(x, on_boundary):
   return on_boundary

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
bc= DirichletBC(V, g, Boundary)

# Solution variable
u = Function(V)
solve(a == L, u, bc)
File("sol.pvd") << u
