"""
Laplace equation with Nitsche method
-Laplace(u) = f
          u = g
Ref: M. Juntunen and R. Stenberg: Nitsche's method for general boundary conditions.
Note: gamma in this code is reciprocal of gamma in the paper
"""

from firedrake import *

mesh = UnitSquareMesh(20,20)

V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

n = FacetNormal(mesh)
h = CellSize(mesh)

u_n = dot(grad(u), n)
v_n = dot(grad(v), n)

gamma = 100.0
f  = Constant(1.0)
g  = Constant(0.0)

a = dot(grad(u),grad(v))*dx - u_n*v*ds - u*v_n*ds + (gamma/h)*u*v*ds
L = f*v*dx + (gamma/h)*g*v*ds - g*v_n*ds

u = Function(V)
solve(a == L, u)
VTKFile("u.pvd").write(u)
