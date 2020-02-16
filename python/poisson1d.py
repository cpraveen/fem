import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from scipy.linalg import solve
import argparse

# See section on "Gauss-Lobattto rules" here
#    https://en.wikipedia.org/wiki/Gaussian_quadrature
def gauss_lobatto_points(n):
    if n == 2:
        x = [-1,1]
    elif n == 3:
        x = [-1,0,1]
    elif n == 4:
        x = [-1, -1/np.sqrt(5), 1/np.sqrt(5), 1]
    elif n == 5:
        x = [-1, -np.sqrt(3.0/7.0), 0, np.sqrt(3.0/7.0), 1]
    elif n == 6:
        x = [-1, -np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
                 -np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 - 2.0*np.sqrt(7)/21), 
                  np.sqrt(1.0/3.0 + 2.0*np.sqrt(7)/21), 
              1]
    elif n == 7:
        x = [-1, -np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
                 -np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  0,
                  np.sqrt(5.0/11.0 - (2.0/11)*np.sqrt(5.0/3.0)),
                  np.sqrt(5.0/11.0 + (2.0/11)*np.sqrt(5.0/3.0)),
              1]
    else:
        print('Not implemented')
        exit()
    return np.array(x)
#------------------------------------------------------------------------------
kw = 10

# rhs function
def f(x):
    return (kw*np.pi)**2 * np.sin(kw*np.pi*x)

# exact solution, use for bc
def uexact(x):
    return x + np.sin(kw*np.pi*x)
#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-degree', type=int, help='Polynomial degree', default=1)
parser.add_argument('-nelem', type=int, help='No. of elements', default=10)
args = parser.parse_args()

xmin, xmax = 0.0, 1.0
N = args.nelem   # number of elements
k = args.degree  # degree of polynomials

xgrid = np.linspace(xmin,xmax,N+1) # grid
hgrid = xgrid[1:] - xgrid[0:-1]    # element lengths

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

# Quadrature points
Nq = k + 1
xq,wq = np.polynomial.legendre.leggauss(Nq)
xq = 0.5*(1 + xq) # transform to [0,1]
wq = 0.5*wq       # transform to [0,1]
print('Number of quadrature points =', Nq)
print('Quadrature points  = ', xq)
print('Quadrature weights = ', wq)

# Evaluate basis function/gradient at quadrature points
shape_value = np.zeros((Nq,k+1))
shape_grad  = np.zeros((Nq,k+1))
for i in range(k+1):
    shape_value[:,i] = shape_funs[i](xq)
    shape_grad [:,i] = shape_grads[i](xq)

b = np.zeros((M,1)) # rhs vector
A = np.zeros((M,M)) # system matrix

# Assemble matrix and rhs
print('Assembling ...')
for n in range(N):
    Aloc = np.zeros((k+1,k+1))
    bloc = np.zeros((k+1,1))
    xphy = xgrid[n] + xq * hgrid[n]
    rhs_values = f(xphy) # rhs function values
    for i in range(k+1):
        bloc[i] = hgrid[n] * np.sum(rhs_values * shape_value[:,i] * wq)
        for j in range(k+1):
            Aloc[i,j] = np.sum(shape_grad[:,i] * shape_grad[:,j] * wq)/hgrid[n]

    # Add to global matrix and vector
    for i in range(k+1):
        ig = local_to_global[n,i]
        b[ig] += bloc[i]
        for j in range(k+1):
            jg = local_to_global[n,j]
            A[ig,jg] += Aloc[i,j]

# Apply bc
print('Appkying bcs ...')
u = np.zeros((M,1))
u[0] = uexact(xgrid[0]); u[-1] = uexact(xgrid[-1])
b -= A@u
b[0] = A[0,0]*u[0]; b[-1] = A[-1,-1]*u[-1]
A[0,1:] = 0.0 # first row
A[1:,0] = 0.0 # first column
A[-1,0:-1] = 0.0  # last row
A[0:-1,-1] = 0.0  # last column

print('Solving Ax=b ...')
u = solve(A,b)

# plot solution
if k == 1:
    xfine = np.linspace(xmin,xmax,1000)
    plt.plot(xfine,uexact(xfine),'k--',xgrid,u,'ro-')
else:
    xfine = np.linspace(xmin,xmax,1000)
    plt.plot(xfine,uexact(xfine),'k--')
    # Sample fem solution on nu uniform points in each element
    nu = 20
    xu = np.linspace(0,1,nu)
    # Create Vandermonde matrix for uniform points
    shape_value = np.zeros((nu,k+1))
    for i in range(k+1):
        shape_value[:,i] = shape_funs[i](xu)
    for n in range(N):
        uloc = u[local_to_global[n,:]] # extract dofs on this element
        values = shape_value @ uloc
        xloc = xgrid[n] + xu * hgrid[n]
        plt.plot(xloc,values,'r-')

plt.legend(('Exact','FEM'))
plt.xlabel('x'); plt.ylabel('u')
plt.title('Degree = '+str(k)+', elements = '+str(N))
plt.show()
