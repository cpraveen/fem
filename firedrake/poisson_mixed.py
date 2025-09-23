# Mixed formulation for the Poisson equation
# See https://www.firedrakeproject.org/demos/poisson_mixed.py.html

from firedrake import *
mesh = UnitSquareMesh(32, 32)

BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)
W = BDM * DG

sigma, u = TrialFunctions(W)
tau, v = TestFunctions(W)

x, y = SpatialCoordinate(mesh)
f = Function(DG).interpolate(
    10*exp(-(pow(x - 0.5, 2) + pow(y - 0.5, 2)) / 0.02))

a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
L = - f*v*dx

bc0 = DirichletBC(W.sub(0), as_vector([0.0, -sin(5*x)]), 3)
bc1 = DirichletBC(W.sub(0), as_vector([0.0, sin(5*x)]), 4)

w = Function(W)

solve(a == L, w, bcs=[bc0, bc1])
sigma, u = w.subfunctions

VTKFile("poisson_mixed.pvd").write(u)

try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

try:
  from firedrake.pyplot import tripcolor
  fig, axes = plt.subplots()
  colors = tripcolor(u, axes=axes)
  fig.colorbar(colors)
  plt.savefig("poisson_mixed.png")
except Exception as e:
  warning("Cannot plot figure. Error msg '%s'" % e)

try:
  plt.show()
except Exception as e:
  warning("Cannot show figure. Error msg '%s'" % e)
