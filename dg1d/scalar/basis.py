from numpy import sqrt

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

def shape_value(n, x):
    return Legendre(n,x)*sqrt(2*n+1)

# Derivatives of Legendre polynomials
def shape_grad(n, x):
    return dLegendre(n,x)*sqrt(2*n+1)
