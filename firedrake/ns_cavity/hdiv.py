# H(div) method for incompressible NS
from firedrake import *

Tf = 1.0
n  = 100     # n x n cells
Re = 100.0   # Reynolds number
Cip = 100.0  # Interior penalty parameter

Cip = Constant(Cip)
mu = Constant(1.0/Re)

# Boundary ids
left, right = 1, 2
bottom, top = 3, 4

mesh = UnitSquareMesh(n,n,"crossed")

BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)
W = BDM * DG

x, y = SpatialCoordinate(mesh)
n  = FacetNormal(mesh)
h  = CellDiameter(mesh)
hf = (h('-') + h('+'))/2.0

# Boundary condition on velocity
vb = conditional(y-1.0 > -1.0e-13, as_vector([1.0,0.0]), as_vector([0.0,0.0]))

def ojump(v, n):
    return outer(v('+'), n('+')) + outer(v('-'), n('-'))

# Parabolic terms
def a(v, vt):
    form =   mu * inner(grad(v), grad(vt))*dx            \
           - mu * inner(ojump(vt,n), avg(grad(v)))*dS  \
           - mu * inner(ojump(v,n), avg(grad(vt)))*dS  \
           + (Cip*mu/hf) * inner(jump(v), jump(vt))*dS \
           - mu * inner(outer(vt,  n), grad(v))*ds       \
           - mu * inner(outer(v-vb,n), grad(vt))*ds      \
           + (Cip*mu/hf) * inner(v-vb, vt)*ds
    return form

# Divergence condition
def b(v, pt):
    form = -div(v)*pt*dx
    return form

# Convective terms
def c(w, v, vt):
    form = -inner(outer(v,w), grad(vt))*dx                \
           + inner(avg(outer(v,w)), ojump(vt,n))*dS            \
           + 0.5*abs(dot(w('+'),n('+')))*inner(jump(v),jump(vt))*dS \
           + dot(vb,vt)*dot(w,n)*ds
    return form

N = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), N])

dt = 0.1
idt = Constant(1.0/dt)

# Define weak formulation
v,  p  = TrialFunctions(W)
vt, pt = TestFunctions(W)
vold = Function(BDM)

F = idt*inner(v - vold, vt)*dx + a(v,vt) + b(vt,p) + c(vold,v,vt) \
    + b(v,pt)

L, f = lhs(F), rhs(F)

# Add time integration loop to reach steady state solution
vtk = VTKFile("sol.pvd")
w = Function(W)
t = 0.0
while t < Tf:
    v, p = w.subfunctions # Unpack deep views that share memory with w
    v.rename("Velocity"); p.rename("Pressure")
    vold.assign(v)
    solve(L == f, w, nullspace=nullspace)
    vtk.write(v)
    t += dt
    print("Time = ", t)