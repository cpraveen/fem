import numpy as np

# f = u^2/2
def flux(x,u):
    return 0.5*u**2

def central_flux(x, ul, ur):
    return 0.5 * (flux(x,ul) + flux(x,ur))

# Godunov flux
def godunov_flux(x, ul, ur):
    u1 = max(0.0, ul)
    u2 = min(0.0, ur)
    return max(flux(x,u1),flux(x,u2))

def roe_flux(x, ul, ur):
    a = np.abs((ul + ur)/2)
    return 0.5 * (flux(x,ul) + flux(x,ur)) - 0.5 * a * (ur - ul)

# Max speed based on cell average values
def max_speed(u):
    return np.abs(u[:,0]).max()
