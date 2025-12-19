"""
This demo program solves the steady incompressible Navier-Stokes equations
for lid-driven cavity problem using Taylor-Hood elements.
Author: Praveen. C
www   : http://math.tifrbng.res.in/~praveen
"""

from firedrake import *
from firedrake.petsc import PETSc

n = 100
Re= 100

# Boundary ids
left, right = 1, 2
bottom, top = 3, 4

# Load mesh from file
mesh = UnitSquareMesh(n,n,"crossed")

# Define function spaces (P2-P1)
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)
W = V * Q

# Define test functions
v,q = TestFunctions(W)

# Define trial functions
w   = Function(W)
u,p = split(w)

# Set parameter values
nu = Constant(1.0/Re)

# Define boundary conditions
bcs = [DirichletBC(W.sub(0), Constant((0,0)), left),
       DirichletBC(W.sub(0), Constant((0,0)), right),
       DirichletBC(W.sub(0), Constant((0,0)), bottom),
       DirichletBC(W.sub(0), Constant((1,0)), top)]

N = VectorSpaceBasis(constant=True, comm=COMM_WORLD)
nullspace = MixedVectorSpaceBasis(W, [W.sub(0), N])

# Weak formulation: F = 0
F =   inner(grad(u)*u, v)*dx \
    + nu*inner(grad(u), grad(v))*dx \
    - div(v)*p*dx \
    - q*div(u)*dx

# Solve F = 0
try:
    solve(F == 0, w, bcs=bcs, nullspace=nullspace,
          solver_parameters={"snes_monitor": None,
                             "ksp_type": "gmres",
                             "mat_type": "aij",
                             "pc_type": "lu",
                             "pc_factor_mat_solver_type": "mumps"})
except PETSc.Error as e:
    if e.ierr == 92:
        warning("MUMPS not installed, skipping direct solve")
    else:
        raise e

u,p = w.subfunctions
u.rename("Velocity"); p.rename("Pressure")

# Compute vorticity by L2 projection
r = TrialFunction(Q)
s = TestFunction(Q)
a = r*s*dx
L = (u[0].dx(1) - u[1].dx(0))*s*dx
vort = Function(Q, name="Vorticity")
solve(a == L, vort)

# Compute stream function
# Laplace(psi) = -vort
a = inner(grad(r), grad(s))*dx
L = vort*s*dx
psi = Function(Q, name="Streamfunction")
wall = DirichletBC(Q, 0, "on_boundary")
solve(a == L, psi, wall)

# Save to file
VTKFile("cavity.pvd").write(u,p,vort,psi)
