"""
-Laplace(u) = 1 in (0,1)x(0,1)
      du/dn = -1/4 on boundary
"""
from dolfin import *

def Boundary(x, on_boundary):
   return on_boundary

n = 50
mesh = UnitSquare(n, n)

V = FunctionSpace(mesh, 'CG', 1)

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
w = Function(V)

solve(a == L, w, bc)

plot(mesh)
plot(w)
interactive()
