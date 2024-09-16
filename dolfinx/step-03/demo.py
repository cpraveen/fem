"""
-Laplace(u) = f in (0,1)x(0,1)
         u  = g on boundary
Choose exact solution as u(x,y) = sin(2*pi*x) cos(2*pi*y) so that the rhs
is given by f(x,y) = (2*pi)^2 u(x,y). The BC g is obtained from exact solution.
"""
from mpi4py import MPI
import dolfinx
from dolfinx import mesh
from dolfinx.fem import functionspace,Function,dirichletbc
from dolfinx.fem import locate_dofs_topological,form,assemble_scalar
from dolfinx.fem.petsc import LinearProblem
from ufl import TrialFunction,TestFunction,dx,inner,grad
from dolfinx.io import VTKFile
from numpy import pi,sin,cos,sqrt

degree = 1
N = [20, 40, 80, 160, 320]

conv = []
for n in N:
    domain = mesh.create_unit_square(MPI.COMM_WORLD, n, n, 
                                    mesh.CellType.triangle)
    V = functionspace(domain, ('CG', degree))
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
                            petsc_options={"ksp_type": "preonly", 
                                           "pc_type": "lu"})
    u = problem.solve()
    u.name = 'u'

    # Interpolate exact solution to higher degree space
    V2 = functionspace(domain, ("CG", degree+1))
    uex = Function(V2)
    uex.interpolate(lambda x: sin(2*pi*x[0])*cos(2*pi*x[1]))

    # Compute error norm
    error_L2 = form(inner(u - uex, u - uex) * dx)
    error_local = assemble_scalar(error_L2)
    error_L2 = sqrt(domain.comm.allreduce(error_local, op=MPI.SUM))

    # Compute hmax
    num_cells = domain.topology.index_map(tdim).size_local
    h = dolfinx.cpp.mesh.h(domain, tdim, range(num_cells))
    hmax_local = h.max()
    hmax = domain.comm.allreduce(hmax_local, op=MPI.MAX)

    error_H1 = error_L2 # TODO: Implement H1 error norm
    print("n = ", n, " h =", hmax, " error = ", error_L2, error_H1)
    conv.append([n, hmax, error_L2, error_H1])

    # Save solution to file
    filename = "output_"+str(n)+".pvd"
    with VTKFile(domain.comm, filename, "w") as vtk:
        vtk.write([u._cpp_object])
    print("Wrote file ", filename)

# Compute convergence rate
from math import log
print("---------------------------------------")
for j in range(5):
   if j==0:
      fmt='{0:4d} {1:14.6e} {2:14.6e} {3:14.6e}'
      print(fmt.format(conv[j][0], conv[j][1], conv[j][2], conv[j][3]))
   else:
      rate_L2 = log(conv[j-1][2]/conv[j][2])/log(2)
      rate_H1 = log(conv[j-1][3]/conv[j][3])/log(2)
      fmt='{0:4d} {1:14.6e} {2:14.6e} {3:14.6e} {4:10.3f} {5:10.3f}'
      print(fmt.format(conv[j][0], conv[j][1], conv[j][2], conv[j][3], rate_L2, rate_H1))
print("---------------------------------------")
