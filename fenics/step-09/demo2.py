"""
-Laplace(u) = 0                 in [-1/2,+1/2] x [0,1]
         u  = (1/pi)*atan2(y,x) on boundary
Experiments to try:
   0) Run test_adapt.py to see the effect of refining a triangle
   1) Set nstep=5 and refine_type='uniform' and run the code. See the solution
      using paraview. Copy conv.dat to conv1.dat
   2) Set nstep=20 and refine_type='adaptive' and run the code. See the solution
      using paraview. Copy conv.dat to conv2.dat
   3) Start matlab and run the matlab code conv.m to see the convergence of L2
      error wrt no. of degrees of freedom.
"""
from dolfin import *
import numpy

# Characteristic function for dirichlet boundary
def Boundary(x, on_boundary):
   return on_boundary

degree = 1

# Exact solution
ue = Expression('(1.0/pi) * atan2(x[1], x[0])',degree=degree+3)

# Initial mesh
n = 10
mesh = RectangleMesh(Point(-0.5, 0), Point(+0.5, 1.0), n, n)

# Number of refinement steps
nstep= 5

# Fraction of cells to refine
REFINE_FRACTION=0.1

# Refinement type: 'uniform' or 'adaptive'
refine_type = 'adaptive'

V = FunctionSpace(mesh, 'CG', degree)

u = TrialFunction(V)
v = TestFunction(V)

# Bilinear form
a = inner(grad(u), grad(v))*dx

# Linear functional
f = Constant(0.0)
L = f*v*dx

# Dirichlet bc
bc = DirichletBC(V, ue, Boundary)

# Solution variable
u = Function(V)

conv = []
file = File('sol.pvd')
ferr = File('error.pvd')
for j in range(nstep):
   solve(a == L, u, bc)
   file << u

   # Compute error in solution
   err = interpolate(ue, V); err.rename("Error","Error")
   err.vector()[:] -= u.vector().array()
   ferr << err
   error_L2 = errornorm(ue, u, norm_type='L2', degree_rise=3)
   error_H1 = errornorm(ue, u, norm_type='H10', degree_rise=3)
   conv.append([V.dim(), mesh.hmax(), error_L2, error_H1])

   n = FacetNormal(mesh)
   dudn = dot( grad(u), n)
   Z = FunctionSpace(mesh, 'DG', 0)
   z = TestFunction(Z)
   ETA = assemble(2*avg(z)*jump(dudn)**2*dS)
   eta = numpy.array([0.5*numpy.sqrt(c.h()*ETA[c.index()]) \
                      for c in cells(mesh)])
   gamma = sorted(eta, reverse=True)[int(len(eta)*REFINE_FRACTION)]
   flag = CellFunction("bool", mesh)
   for c in cells(mesh):
      flag[c] = eta[c.index()] > gamma

   # refine the mesh
   if j < nstep-1:
      if refine_type == 'adaptive':
         mesh_new = adapt(mesh, flag)
      else:
         mesh_new = adapt(mesh)

   V = adapt(V, mesh_new)
   u = adapt(u, mesh_new)
   a = adapt(Form(a), mesh_new)
   L = adapt(Form(L), mesh_new)
   bc= adapt(bc, mesh_new)
   mesh = mesh_new



print "---------------------------------------"
f = open('conv.dat','w')
for j in range(nstep):
   fmt='{0:6d} {1:14.6e} {2:14.6e} {3:14.6e}'
   print fmt.format(conv[j][0], conv[j][1], conv[j][2], conv[j][3])
   f.write(str(conv[j][0])+' '+str(conv[j][1])+' '+str(conv[j][2]))
   f.write(' '+str(conv[j][3])+'\n')
print "---------------------------------------"
f.close()
