"""
-Laplace(u) = 1 in (0,1)x(0,1)
      du/dn = -1/4 on boundary
"""
from dolfin import *

n = 50
mesh = UnitSquareMesh(n, n)
degree = 1
V = FunctionSpace(mesh, 'CG', degree)
u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Constant(1)
g = Constant(-0.25)
L = f*v*dx + g*v*ds

bc= DirichletBC(V, 0, 'near(x[0],0) && near(x[1],0)', 'pointwise')

# Solution variable
u = Function(V)
solve(a == L, u, bc)
u.rename('solution','solution')
File("sol.pvd") << u