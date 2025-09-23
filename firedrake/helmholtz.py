# Simple Helmholtz equation
# See https://www.firedrakeproject.org/demos/helmholtz.py.html

from firedrake import *
mesh = UnitSquareMesh(10, 10)

V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

f = Function(V)
x, y = SpatialCoordinate(mesh)
f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))

a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
L = inner(f, v) * dx

u = Function(V)
solve(a == L, u, solver_parameters={'ksp_type': 'cg', 'pc_type': 'none'})
VTKFile("helmholtz.pvd").write(u)

try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

try:
  from firedrake.pyplot import tripcolor, tricontour
  fig, axes = plt.subplots()
  colors = tripcolor(u, axes=axes)
  fig.colorbar(colors)
except Exception as e:
  warning("Cannot plot figure. Error msg: '%s'" % e)

try:
  fig, axes = plt.subplots()
  contours = tricontour(u, axes=axes)
  fig.colorbar(contours)
  plt.savefig("helmholtz.png")
except Exception as e:
  warning("Cannot plot figure. Error msg: '%s'" % e)

try:
  plt.show()
except Exception as e:
  warning("Cannot show figure. Error msg: '%s'" % e)

f.interpolate(cos(x*pi*2)*cos(y*pi*2))
print(sqrt(assemble(dot(u - f, u - f) * dx)))
