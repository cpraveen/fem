from numpy import sqrt

# Legendre polynomials on [-1,+1]
def shape_value(i, x):
    if i==0:
        value = 1.0
    elif i==1:
        value = x
    elif i==2:
        value = 0.5*(3*x**2 - 1)
    elif i==3:
        value = 0.5*(5*x**3 - 3*x)
    else:
        print("shape_value not implemented")
        exit()
    return value*sqrt(2*i+1)

# Derivatives of Legendre polynomials
def shape_grad(i, x):
    if i==0:
        value = 0.0
    elif i==1:
        value = 1.0
    elif i==2:
        value = 3*x
    elif i==3:
        value = 0.5*(15*x**2 - 3)
    else:
        print("shape_grad not implemented")
        exit()
    return value*sqrt(2*i+1)
