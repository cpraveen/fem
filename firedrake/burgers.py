# Burgers equation
# See https://www.firedrakeproject.org/demos/burgers.py.html

from firedrake import *
n = 30
mesh = UnitSquareMesh(n, n)

V = VectorFunctionSpace(mesh, "CG", 2)
V_out = VectorFunctionSpace(mesh, "CG", 1)

u_ = Function(V, name="Velocity")
u = Function(V, name="VelocityNext")

v = TestFunction(V)

x = SpatialCoordinate(mesh)
ic = project(as_vector([sin(pi*x[0]), 0]), V)

u_.assign(ic)
u.assign(ic)

nu = 0.0001

timestep = 1.0/n

F = (inner((u - u_)/timestep, v)
     + inner(dot(u,nabla_grad(u)), v) + nu*inner(grad(u), grad(v)))*dx

outfile = VTKFile("burgers.pvd")

outfile.write(project(u, V_out, name="Velocity"))

t = 0.0
end = 0.5
while (t <= end):
    solve(F == 0, u)
    u_.assign(u)
    t += timestep
    outfile.write(project(u, V_out, name="Velocity"))
