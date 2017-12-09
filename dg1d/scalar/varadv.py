import numpy as np

def speed(x):
    return 1.0 + (1.0 - x**2)**5

# f = a(x)*u
def flux(x,u):
    return speed(x)*u

# Upwind flux
def numflux(x, ul, ur):
    a = speed(x)
    return max(a,0.0)*ul + min(a,0.0)*ur

def max_speed(u):
    x = np.linspace(-1.0,1.0,100)
    return np.abs(speed(x)).max()
