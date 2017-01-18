"""
-div( mu * grad(u) ) = 1 in [-1,+1] x [-1,+1]
                   u = 0 on boundary

mu = mu1 for x^2 + y^2 < 0.5^2
     mu2 otherwise
"""
from dolfin import *
import numpy

# Characteristic function for dirichlet boundary
def Boundary(x, on_boundary):
   return on_boundary

mu1 = 1.0
mu2 = 10.0

class Coefficient(Expression):
   def eval(self, values, x):
      if x[0]**2 + x[1]**2 <= 0.5**2:
         values[0] = mu1
      else:
         values[0] = mu2

# Initial mesh
n = 10
mesh = Rectangle(-1.0, -1.0, +1.0, +1.0, n, n)

# Number of refinement steps
nstep= 20

# Fraction of cells to refine
REFINE_FRACTION=0.1

# Refinement type: 'uniform' or 'adaptive'
refine_type = 'adaptive'

file = File('sol.pvd')
for j in range(nstep):
   V = FunctionSpace(mesh, 'CG', 1)

   u = TrialFunction(V)
   v = TestFunction(V)

   # Bilinear form
   mu = Coefficient()
   a = mu*inner(grad(u), grad(v))*dx

   # Linear functional
   f = Constant(1.0)
   L = f*v*dx

   # Dirichlet bc
   bc = DirichletBC(V, 0.0, Boundary)

   # Solution variable
   u = Function(V)

   solve(a == L, u, bc)

   file << u

   n = FacetNormal(mesh)
   dudn = dot( grad(u), n)
   Z = FunctionSpace(mesh, 'DG', 0)
   z = TestFunction(Z)
   ETA = assemble(2*avg(z)*jump(dudn)**2*dS)
   eta = numpy.array([0.5*numpy.sqrt(c.diameter()*ETA[c.index()]) \
                      for c in cells(mesh)])
   gamma = sorted(eta, reverse=True)[int(len(eta)*REFINE_FRACTION)]
   flag = CellFunction("bool", mesh)
   for c in cells(mesh):
      flag[c] = eta[c.index()] > gamma

   # refine the mesh
   if j < nstep-1:
      if refine_type == 'adaptive':
         mesh = refine(mesh, flag)
      else:
         mesh = refine(mesh)
