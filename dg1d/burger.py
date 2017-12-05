# f = u^2/2
def flux(u):
    return 0.5*u**2

# Godunov flux
def numflux(ul, ur):
    u1 = max(0.0, ul)
    u2 = min(0.0, ur)
    return max(flux(u1),flux(u2))
