"""
-Laplace(u) = 0
         u  = (1/pi)*atan2(y,x)
"""
from dolfin import *
from math import atan2, log

def BottomLeft(x, on_boundary):
   return near(x[1],0) and x[0] < DOLFIN_EPS and on_boundary

def BottomRight(x, on_boundary):
   return near(x[1],0) and x[0] > DOLFIN_EPS and on_boundary

# Remaining part of boundary
def Boundary(x, on_boundary):
   return x[1] > DOLFIN_EPS and on_boundary

# Exact solution
ue = Expression('(1.0/pi)*atan2(x[1], x[0])')

mesh = Mesh('mesh.xml')

nstep= 5
conv = []
file = File('sol.pvd')
for j in range(nstep):
   V = FunctionSpace(mesh, 'CG', 1)

   u = TrialFunction(V)
   v = TestFunction(V)

   # Bilinear form
   a = inner(grad(u), grad(v))*dx

   # Linear functional
   f = Constant(0.0)
   L = f*v*dx

   # Dirichlet bc
   bc1 = DirichletBC(V, ue, Boundary)
   bc2 = DirichletBC(V, Constant(1), BottomLeft)
   bc3 = DirichletBC(V, Constant(0), BottomRight)
   bc  = [bc1, bc2, bc3]

   # Solution variable
   w = Function(V)

   solve(a == L, w, bc)

   w.rename("sol","sol")
   file << w
   error_L2 = errornorm(ue, w, norm_type='L2', degree_rise=3)
   conv.append([mesh.hmax(), error_L2])

   # refine the mesh
   mesh = refine(mesh)

# Compute convergence rate
print "---------------------------------------"
for j in range(nstep):
   if j==0:
      fmt='{0:14.6e} {1:14.6e}'
      print fmt.format(conv[j][0], conv[j][1])
   else:
      rate_L2 = log(conv[j-1][1]/conv[j][1])/log(2)
      fmt='{0:14.6e} {1:14.6e} {2:10.3f}'
      print fmt.format(conv[j][0], conv[j][1], rate_L2)
print "---------------------------------------"
