'''
Solve
   -Laplace(u) = f  in  (0,1) x (0,1)
            u  = g  on  boundary
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

mesh = meshio.read("mesh_tri.msh")
x, y = mesh.points[:,0], mesh.points[:,1]
cells = mesh.cells_dict["triangle"]
faces = mesh.cells_dict["line"]

# Find unique boundary points
bpts = np.unique(faces)

npoints = len(x)
ncells = cells.shape[0]
nfaces = faces.shape[0]
nbpts = len(bpts)

print('points, cells, faces, bpts = ', npoints, ncells, nfaces, nbpts)

# RHS function
f = lambda x,y: np.sin(2.0*np.pi*x) * np.cos(2.0*np.pi*y)

# Exact solution
uexact = lambda x,y: (0.5/(2.0*np.pi)**2) * np.sin(2.0*np.pi*x) * np.cos(2.0*np.pi*y)

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

# Fill boundary values
for i in range(nbpts):
    p = bpts[i]
    u[p] = uexact(x[p], y[p])

# Modify matrix and rhs to apply bc
b -= A @ u
for i in range(nbpts):
    c = bpts[i]
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

plt.show()
