"""
-Laplace(u) = f in (0,1)x(0,1)
         u  = g on boundary
Choose exact solution as u(x,y) = sin(2*pi*x) cos(2*pi*y) so that the rhs 
is given by f(x,y) = 2(2*pi)^2 u(x,y). The BC g is obtained from exact solution.
"""
from dolfin import *

degree = 1
N = [20, 40, 80, 160, 320]

conv = []
file = File('sol.pvd')
for n in N:
   mesh = UnitSquareMesh(n,n)

   V = FunctionSpace(mesh, 'CG', degree)
   u = TrialFunction(V)
   v = TestFunction(V)

   # Bilinear form
   a = inner(grad(u), grad(v))*dx

   # Linear functional
   f = Expression('8*pi*pi*sin(2*pi*x[0])*cos(2*pi*x[1])',degree=degree)
   L = f*v*dx

   # Dirichlet bc
   g = Expression('sin(2*pi*x[0])*cos(2*pi*x[1])',degree=degree)
   bc= DirichletBC(V, g, 'on_boundary')

   # Solution variable
   u = Function(V)
   solve(a == L, u, bc)
   u.rename('u','u')
   file << u

   uexact = Expression('sin(2*pi*x[0])*cos(2*pi*x[1])',degree=degree+3)
   error_L2 = errornorm(uexact, u, norm_type='L2')
   error_H1 = errornorm(uexact, u, norm_type='H10')
   print("n = ", n, " h =", mesh.hmax(), " error = ", error_L2, error_H1)
   conv.append([n, mesh.hmax(), error_L2, error_H1])

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
