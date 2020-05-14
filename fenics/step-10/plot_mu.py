from dolfin import *
import matplotlib.pyplot as plt

mu1 = 1.0
mu2 = 10.0

class Coefficient(UserExpression):
   def eval(self, values, x):
      if x[0]**2 + x[1]**2 <= 0.5**2:
         values[0] = mu1
      else:
         values[0] = mu2
   def value_shape(self):
         return ()

n = 10
mesh = RectangleMesh(Point(-1.0, -1.0), Point(+1.0, +1.0), n, n)
mu = Coefficient(degree=0)
V0 = FunctionSpace(mesh,'DG',0)
mu = project(mu,V0)
c = plot(mu)
plt.colorbar(c)
plt.savefig('mu.pdf')
print('Coefficient mu printed to mu.pdf')
