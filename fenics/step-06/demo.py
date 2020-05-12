"""
-Laplace(u) = 0
         u  = (1/pi)*atan2(y,x)
"""
from dolfin import *
from math import atan2, log

# x < 0 and y = 0
def BottomLeft(x, on_boundary):
   return near(x[1],0) and x[0] <= 0.0 and on_boundary

# x > 0 and y = 0
def BottomRight(x, on_boundary):
   return near(x[1],0) and x[0] > 0.0 and on_boundary

# Remaining part of boundary
def Boundary(x, on_boundary):
   return x[1] > 0.0 and on_boundary

# Polynomial degree
degree = 1

# Exact solution
ue = Expression('(1.0/pi)*atan2(x[1], x[0])',degree=degree+3)

mesh = Mesh('mesh.xml')

nstep= 5
conv = []
file = File('sol.pvd')
for j in range(nstep):
   V = FunctionSpace(mesh, 'CG', degree)

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

   # We can use exact solution for all dirichlet bc. Comment above line
   # and uncomment following line.
   # bc  = DirichletBC(V, ue, 'on_boundary')

   # Solution variable
   u = Function(V)
   solve(a == L, u, bc)

   u.rename("sol","sol")
   file << u
   error_L2 = errornorm(ue, u, norm_type='L2')
   error_H1 = errornorm(ue, u, norm_type='H10')
   conv.append([mesh.hmax(), error_L2, error_H1])

   # refine the mesh
   mesh = refine(mesh)

# Compute convergence rate
print("----------------------------------------------------------------------")
fmt='{0:14.6e} {1:14.6e} {2:10.3f} {3:14.6e} {4:10.3f}'
for j in range(nstep):
   if j==0:
      print(fmt.format(conv[j][0], conv[j][1], 0, conv[j][2], 0))
   else:
      rate_L2 = log(conv[j-1][1]/conv[j][1])/log(2)
      rate_H1 = log(conv[j-1][2]/conv[j][2])/log(2)
      print(fmt.format(conv[j][0], conv[j][1], rate_L2, conv[j][2], rate_H1))
print("----------------------------------------------------------------------")
