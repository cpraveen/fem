"""
-Laplace(u) = 1 in (0,1)x(0,1)
         u  = 0 on boundary
"""
from mpi4py import MPI
from dolfinx import mesh
from dolfinx.fem import FunctionSpace,Function,Constant,dirichletbc
from dolfinx.fem import locate_dofs_topological
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction,TestFunction,dx,inner,grad
from petsc4py.PETSc import ScalarType
from dolfinx.io import VTKFile,XDMFFile

domain = mesh.create_unit_square(MPI.COMM_WORLD, 20, 20, 
                                 mesh.CellType.triangle)
degree = 1
V = FunctionSpace(domain, ('CG', degree))
u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Constant(domain, ScalarType(1))
L = f*v*dx

# Create facet to cell connectivity required to determine boundary facets
tdim = domain.topology.dim
fdim = tdim - 1
domain.topology.create_connectivity(fdim, tdim)
boundary_facets = mesh.exterior_facet_indices(domain.topology)

# Dirichlet bc
g = Constant(domain, ScalarType(0))
boundary_dofs = locate_dofs_topological(V, fdim, boundary_facets)
bc = dirichletbc(g, boundary_dofs, V)

# Define linear problem and solve
problem = LinearProblem(a, L, bcs=[bc], 
                        petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
u = problem.solve()
u.name = 'u'

# Save solution to file
with VTKFile(domain.comm, "output.pvd", "w") as vtk:
    vtk.write([u._cpp_object])
with XDMFFile(domain.comm, "output.xdmf", "w") as xdmf:
    xdmf.write_mesh(domain)
    xdmf.write_function(u)
