'''
Solve
   -Laplace(u)  = 1,   {1 < r < 2}
            u   = 0,   {r = 1}
          du/dn = 0,   {r = 2}
'''
import meshio
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
import numpy.linalg as nla

# Assemble matrix and rhs on one triangle
def assemble_cell(x, y, f):
    o = np.ones(3)
    H = np.array([o, x, y])
    area = 0.5 * nla.det(H)
    G = nla.inv(H) @ np.array([[0.0, 0.0],
                               [1.0, 0.0],
                               [0.0, 1.0]])
    A = 0.5 * area * G @ G.T
    xc = np.sum(x) / 3.0
    yc = np.sum(y) / 3.0
    b = np.zeros(3)
    b[:] = (1.0/6.0) * area * f(xc, yc)
    return A, b

mesh = meshio.read("annulus.msh")
x, y = mesh.points[:,0], mesh.points[:,1]
cells = mesh.cells_dict["triangle"]
allfaces = mesh.cells_dict["line"]
ifaces = mesh.cell_sets_dict["inner"]["line"]
faces = allfaces[ifaces]

# Find unique dirichlet boundary points
bpts = np.unique(faces)

npoints = len(x)
ncells = cells.shape[0]
nfaces = faces.shape[0]
nbpts = len(bpts)

print('points, cells, faces, bpts = ', npoints, ncells, nfaces, nbpts)

# RHS function
f = lambda x,y: 1.0

# boundary value on inner boundary
ubdry = lambda x,y: 0.0

b = np.zeros(npoints)   # rhs vector
A = lil_matrix((npoints,npoints)) # system matrix

# Loop over cells and assemble
for c in range(ncells):
    pts = cells[c,:]
    xv, yv = x[pts], y[pts]
    Aloc, bloc = assemble_cell(xv, yv, f)
    # Copy local to global
    for i in range(3):
        ig = cells[c,i]
        b[ig] += bloc[i]
        for j in range(3):
            jg = cells[c,j]
            A[ig,jg] += Aloc[i,j]

# Solution array
u = np.zeros(npoints)

# Fill boundary values into solution array u
for p in bpts:
    u[p] = ubdry(x[p], y[p])

# Modify matrix and rhs to apply dirichlet bc
b -= A @ u
for c in bpts:
    b[c] = A[c,c] * u[c]
    A[c, :c]     = 0.0 # other entries in c'th row
    A[c, (c+1):] = 0.0
    A[:c, c]     = 0.0 # other entries in c'th column
    A[(c+1):,c]  = 0.0

A = csc_matrix(A) # convert to csc since spsolve needs this

print('Solving Ax=b ...')
u = spsolve(A, b)

# Plot solution
plt.figure()
c = plt.tricontourf(x, y, cells, u, cmap='rainbow', levels=20)
plt.colorbar(c)
plt.tricontour(x, y, cells, u, colors='k', levels=20)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution')
plt.axis('equal')
plt.savefig('sol.svg')

uexact = lambda r: -0.25 * r**2 + 2 * np.log(r) + 0.25

j1 = np.where(np.abs(y) < 1.0e-10)
j2 = np.where(x > 0)
j  = np.intersect1d(j1,j2)
plt.figure()
plt.plot(x[j], u[j], 'o', label='FEM')
x1 = np.sort(x[j])
plt.plot(x1, uexact(x1), '-', label='Exact')
plt.xlabel('r')
plt.ylabel('u(r)')
plt.legend()
plt.title('Solution along 1 <= x <= 2, y = 0')
plt.savefig('line.svg')

plt.show()
