from numpy import sqrt,empty

# Legendre polynomials on [-1,+1]
def Legendre(n, x):
    if n==0:
        value = 1.0
    elif n==1:
        value = x
    else:
        value = (2.0*n-1.0)/n * x * Legendre(n-1,x) - (n-1.0)/n * Legendre(n-2,x)
    return value

# Derivative of Legendre
def dLegendre(n, x):
    if n==0:
        value = 0.0
    elif n==1:
        value = 1.0
    else:
        value = n * Legendre(n-1,x) + x * dLegendre(n-1,x)
    return value

# Returns n'th basis function evaluated at x in [-1,+1]
def shape_value(n, x):
    return Legendre(n,x)*sqrt(2*n+1)

# Returns derivative of n'th basis function evaluated at x in [-1,+1]
def shape_grad(n, x):
    return dLegendre(n,x)*sqrt(2*n+1)

# Vandermonde matrix for degree = k at points defined in x
# V_ij = phi_j( x_i )
def Vandermonde(k, x):
    n = len(x)
    V = empty((n, k+1))
    for i in range(n):
        for j in range(k+1):
            V[i, j] = shape_value(j, x[i])
    return V
