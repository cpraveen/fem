"""
-Laplace(u) = f in (0,1)x(0,1)
         u  = g on boundary
Choose exact solution as u(x,y) = sin(2*pi*x) cos(2*pi*y) so that the rhs
is given by f(x,y) = (2*pi)^2 u(x,y). The BC g is obtained from exact solution.
"""
from mpi4py import MPI
from dolfinx import mesh
from dolfinx.fem import FunctionSpace,Function,Constant,dirichletbc
from dolfinx.fem import locate_dofs_topological
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction,TestFunction,dx,inner,grad,SpatialCoordinate
from petsc4py.PETSc import ScalarType
from dolfinx.io import VTKFile,XDMFFile
from numpy import pi,sin,cos

domain = mesh.create_unit_square(MPI.COMM_WORLD, 20, 20, 
                                 mesh.CellType.triangle)
degree = 1
V = FunctionSpace(domain, ('CG', degree))
u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Function(V)
f.interpolate(lambda x: 8*pi*pi*sin(2*pi*x[0])*cos(2*pi*x[1]))
L = f*v*dx

# Create facet to cell connectivity required to determine boundary facets
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)

# Dirichlet bc
g = Function(V)
g.interpolate(lambda x: sin(2*pi*x[0])*cos(2*pi*x[1]))
boundary_dofs = locate_dofs_topological(V, fdim, boundary_facets)
bc = dirichletbc(g, boundary_dofs)

# Define linear problem and solve
problem = LinearProblem(a, L, bcs=[bc], 
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
u = problem.solve()
u.name = 'u'

# Save solution to file
with VTKFile(domain.comm, "output.pvd", "w") as vtk:
    vtk.write([u._cpp_object])
