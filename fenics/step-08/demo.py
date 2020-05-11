"""
-Laplace(u) = f in (0,1)x(0,1)
         u  = g on boundary
Choose exact solution as u(x,y) = sin(2*pi*x) cos(2*pi*y) so that the rhs 
is given by f(x,y) = (2*pi)^2 u(x,y). The BC g is obtained from exact solution.
"""
from dolfin import *
import numpy

degree = 1
mesh = UnitSquareMesh(20,20)

V = FunctionSpace(mesh, 'CG', degree)

u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx
A = assemble(a)

# Linear functional
f = Expression('8*pi*pi*sin(2*pi*x[0])*cos(2*pi*x[1])',degree=degree+3)
L = f*v*dx
b = assemble(L)

# Dirichlet bc
g = Expression('sin(2*pi*x[0])*cos(2*pi*x[1])',degree=degree+3)
bc= DirichletBC(V, g, 'on_boundary')
bc.apply(A, b)

# Solution variable
u = Function(V)
solve(A, u.vector(), b)
File('u.pvd') << u

# project gradient
V2 = VectorFunctionSpace(mesh, 'CG', degree)
w = TrialFunction(V2)
v = TestFunction(V2)
ag = inner(w, v)*dx
Lg = inner(grad(u), v)*dx
gradu = Function(V2)
solve(ag == Lg, gradu)

# Alternately, we can use the "project" function
#gradu = project(grad(u), V2)

File('gradu.pvd') << gradu

# Maximum norm computation
ue = Function(V)
ue = interpolate(g, V)
# Get numpy arrays
u_array = u.vector().get_local()
ue_array= ue.vector().get_local()
print("Max error =", numpy.abs(u_array - ue_array).max())
