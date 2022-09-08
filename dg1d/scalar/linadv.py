# f = a*u with a = 1
def flux(x,u):
    return u

# Central flux
def central_flux(x, ul, ur):
    return 0.5 * (ul + ur)

# Upwind flux
def upwind_flux(x, ul, ur):
    return ul

def max_speed(u):
    return 1.0
