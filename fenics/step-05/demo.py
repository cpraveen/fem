"""
-Laplace(u) = 0
         u  = r^(2/3) * sin(2*theta/3) on boundary
"""
from dolfin import *
from math import atan2, sqrt, sin, log

def Boundary(x, on_boundary):
   return on_boundary

# Exact solution
class uexact(Expression):
   def eval(self, values, x):
      r = sqrt(x[0]**2 + x[1]**2)
      theta = atan2(x[1], x[0])
      if theta < 0:
         theta = theta + 2*pi
      values[0] = r**(2.0/3.0) * sin(2.0*theta/3.0)

mesh = Mesh("Gamma.xml")

conv = []
file = File('sol.pvd')
for j in range(5):
   V = FunctionSpace(mesh, 'CG', 1)

   u = TrialFunction(V)
   v = TestFunction(V)

   # Bilinear form
   a = inner(grad(u), grad(v))*dx

   # Linear functional
   f = Constant(0.0)
   L = f*v*dx

   # Dirichlet bc
   ue = uexact()
   bc = DirichletBC(V, ue, Boundary)

   # Solution variable
   w = Function(V)

   solve(a == L, w, bc)

   file << w
   error_L2 = errornorm(ue, w, norm_type='L2', degree_rise=3)
   conv.append([V.dim(), mesh.hmax(), error_L2])

   # refine the mesh
   mesh = refine(mesh)


print "---------------------------------------"
f = open('conv.dat','w')
for j in range(5):
   if j==0:
      fmt='{0:6d} {1:14.6e} {2:14.6e}'
      print fmt.format(conv[j][0], conv[j][1], conv[j][2])
   else:
      rate_L2 = log(conv[j-1][2]/conv[j][2])/log(2)
      fmt='{0:6d} {1:14.6e} {2:14.6e} {3:10.3f}'
      print fmt.format(conv[j][0], conv[j][1], conv[j][2], rate_L2)
   f.write(str(conv[j][0])+' '+str(conv[j][1])+' '+str(conv[j][2])+'\n')
print "---------------------------------------"
f.close()
