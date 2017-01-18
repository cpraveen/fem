"""
Scalar advection on circles using DG
a.grad(u) = 0 in [0,1]x[0,1],  a = (y, -x)
u = uin on x=0

Since div(a)=0, we write DG scheme for 
div(a*u) = 0
"""
from dolfin import *

class uinlet(Expression):
   def eval(self, value, x):
      if abs(x[1]-0.5) < 0.25:
         value[0]=1.0
      else:
         value[0]=0.0

np = 50

mesh = UnitSquare(np, np)

X = FunctionSpace(mesh, "DG", 1)

v = TestFunction(X)
u = TrialFunction(X)

a = Expression(("x[1]", "-x[0]"))

n   = FacetNormal(mesh)
an  = dot(a,n)
anp = 0.5*(an+abs(an))
anm = 0.5*(an-abs(an))
H   = anp('+')*u('+') + anm('+')*u('-')

form = -inner(a, grad(v))*u*dx + H*jump(v)*dS + v*anp*u*ds + anm*v*uinlet()*ds
a    = lhs(form)
L    = rhs(form)

u = Function(X)
solve(a == L, u, bcs=None)
File("sol.pvd") << u
plot(u, interactive=True)
