'''
Solves the bvp
    -u'' = f in (xmin,xmax)
     u(xmin) = alpha
     u(xmax) = beta
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from scipy.sparse import lil_matrix, csc_matrix
from scipy.sparse.linalg import spsolve
import argparse
from lgl import *

#------------------------------------------------------------------------------
kw = 20

# rhs function
def f(x):
    return (kw)**2 * np.sin(kw*x)

# exact solution, use for bc
def exact_solution(x):
    return x + np.sin(kw*x)

# Returns boundary condition
def boundary_value(x):
    return exact_solution(x)
#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-nelem', type=int, help='No. of elements', default=10)
parser.add_argument('-plot_grid', type=int,
                    help='No. of points for plotting solution', default=20)
args = parser.parse_args()

xmin, xmax = 0.0, 1.0
N = args.nelem   # number of elements
k = args.degree  # degree of polynomials

# Make grid
xgrid = np.linspace(xmin,xmax,N+1) # grid
hgrid = xgrid[1:] - xgrid[0:-1]    # element lengths

# Build local dof to global dof map
local_to_global = np.zeros((N,k+1),dtype=int)
count = 0
for i in range(N): # element loop
    for j in range(k+1): # dofs inside element
        if i>0 and j==0:
            local_to_global[i,j] = local_to_global[i-1,k]
        else:
            local_to_global[i,j] = count
            count += 1

M = count # total unknowns, including boundary values

print('Number of elements = ', N)
print('Degree             = ', k)
print('Number of unknowns = ', M)

# Generate Lagrange polynomials and their derivatives
# WARNING: Dont use for very high degree, < 10 should be ok
xs = 0.5*(1 + gauss_lobatto_points(k+1))
shape_funs, shape_grads = [], []
for i in range(k+1):
    values = np.zeros((k+1,1))
    values[i] = 1.0
    shape_fun = lagrange(xs, values)
    shape_funs.append(shape_fun)
    shape_grads.append(shape_fun.deriv())

# Gauss-Legendre quadrature points
Nq = k + 1
xq,wq = np.polynomial.legendre.leggauss(Nq)
xq = 0.5*(1 + xq) # transform to [0,1]
wq = 0.5*wq       # transform to [0,1]
print('Number of quadrature points =', Nq)

# Evaluate basis function/gradient at quadrature points
shape_value = np.zeros((Nq,k+1))
shape_grad  = np.zeros((Nq,k+1))
for i in range(k+1):
    shape_value[:,i] = shape_funs[i](xq)
    shape_grad [:,i] = shape_grads[i](xq)

b = np.zeros((M,1))   # rhs vector
A = lil_matrix((M,M)) # system matrix

# Assemble matrix and rhs
print('Assembling ...')
for n in range(N): # Loop over elements
    Aloc = np.zeros((k+1,k+1))      # local matrix
    bloc = np.zeros((k+1,1))        # local rhs vector
    xphy = xgrid[n] + xq * hgrid[n] # quad pts in real space
    rhs_values = f(xphy)            # rhs function values
    for i in range(k+1): # Loop over basis functions
        bloc[i] = hgrid[n] * np.sum(rhs_values * shape_value[:,i] * wq)
        for j in range(k+1): # Loop over basis functions
            Aloc[i,j] = np.sum(shape_grad[:,i] * shape_grad[:,j] * wq)/hgrid[n]

    # Add to global matrix and vector
    for i in range(k+1):
        ig = local_to_global[n,i]
        b[ig] += bloc[i]
        for j in range(k+1):
            jg = local_to_global[n,j]
            A[ig,jg] += Aloc[i,j]

# Apply bc
print('Applying bcs ...')
u = np.zeros((M,1))
u[0]  = boundary_value(xgrid[0])
u[-1] = boundary_value(xgrid[-1])
b    -= A@u
b[0]  = A[0,0]*u[0]
b[-1] = A[-1,-1]*u[-1]
A[0,1:]    = 0.0 # first row
A[1:,0]    = 0.0 # first column
A[-1,0:-1] = 0.0 # last row
A[0:-1,-1] = 0.0 # last column
A = csc_matrix(A) # convert to csc since spsolve needs this

print('Solving Ax=b ...')
u = spsolve(A,b)

# Compute error norm

# plot solution
if k == 1:
    xfine = np.linspace(xmin,xmax,1000)
    plt.plot(xfine,exact_solution(xfine),'k--',xgrid,u,'ro-')
else: # sub-sample inside each element
    xfine = np.linspace(xmin,xmax,1000)
    plt.plot(xfine,exact_solution(xfine),'k--')
    # Sample fem solution on nu uniform points in each element
    nu = args.plot_grid
    xu = np.linspace(0,1,nu)
    # Create Vandermonde matrix for uniform points
    shape_value = np.zeros((nu,k+1))
    for i in range(k+1):
        shape_value[:,i] = shape_funs[i](xu)
    for n in range(N): # Loop over elements
        uloc = u[local_to_global[n,:]] # extract dofs on this element
        values = shape_value @ uloc
        xloc = xgrid[n] + xu * hgrid[n]
        plt.plot(xloc,values,'r-')

plt.plot(xgrid,0*xgrid,'|-',markersize=15)
plt.legend(('Exact','FEM'))
plt.xlabel('x'); plt.ylabel('u')
plt.title('Degree = '+str(k)+', elements = '+str(N))
plt.grid(True)
plt.show()
