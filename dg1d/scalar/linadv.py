# f = a*u with a = 1
def flux(x,u):
    return u

# Upwind flux
def numflux(x, ul, ur):
    return ul

def max_speed(u):
    return 1.0
