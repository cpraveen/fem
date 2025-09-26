"""
Solve -Laplace(u) = f in [0,1]x[0,1]
with u = g on boundary
"""
from firedrake import *

k    = 1    # degree of polynomial
Cip  = 10.0 # IP penalty
np   = 50   # number of points on each side of square

mesh = UnitSquareMesh(np,np)
V    = FunctionSpace(mesh, 'DG', k)
u    = TrialFunction(V)
v    = TestFunction(V)
n    = FacetNormal(mesh)
h    = CellDiameter(mesh)
havg = (h('-')+h('+'))/2.0

f = Constant(1.0) # RHS function
g = Constant(0.0) # Dirichlet bc

F =   inner(grad(u),grad(v))*dx         \
    - inner(avg(grad(u)), jump(v,n))*dS \
    - inner(avg(grad(v)), jump(u,n))*dS \
    - inner(grad(u), v*n)*ds            \
    - inner(grad(v), (u-g)*n)*ds        \
    + (Cip*k*k/havg)*jump(u)*jump(v)*dS \
    + (Cip*k*k/h)*(u-g)*v*ds            \
    - f*v*dx

a,L = lhs(F), rhs(F)
u = Function(V)
solve(a==L, u)
VTKFile('sol.pvd').write(u)
