import numpy as np

gasGam = 1.4
gasR   = 1.0

# Compute pressure given density, momentum and energy
def pressure(rho,mom,ene):
    return (gasGam-1.0)*(ene - 0.5*mom**2/rho)

# Compute Euler flux given density, momentum, energy
def flux(rho,mom,ene):
    pre  = pressure(rho,mom,ene)
    frho = mom
    fmom = pre + mom**2/rho
    fene = (ene + pre)*mom/rho
    return frho,fmom,fene

# Compute maximum wave speed
# u = conserved variables
def max_eig(u):
    vel = u[1]/u[0]
    pre = pressure(u[0],u[1],u[2])
    lam = np.abs(vel) + np.sqrt(gasGam * pre / u[0])
    return lam

# Compute maximum wave speed in each cell based on cell average value
def max_speed(rho,mom,ene):
    vel = mom[:,0] / rho[:,0]
    pre = pressure(rho[:,0], mom[:,0], ene[:,0])
    lam = np.abs(vel) + np.sqrt(gasGam * pre / rho[:,0])
    return lam.max()

# Compute matrix of right and left eigenvectors from u = conserved variables
# R = matrix of right eigenvectors
# L = matrix of left eigenvectors
# L = R^{-1} so that R*L = I
def EigMatrix(u):
    g1 = gasGam - 1.0
    g2 = 0.5*g1

    d = u[0]
    v = u[1] / d
    p = pressure(u[0],u[1],u[2])
    c = np.sqrt(gasGam * p / d)
    h = c**2 / g1 + 0.5 * v**2
    f = 0.5 * d / c 
   
    R, L = np.zeros((3,3)), np.zeros((3,3))

    # Inverse eigenvector-matrix
    L[0,0] =  1.0 - g2 * v**2 / c**2
    L[1,0] =  (g2 * v**2 - v * c) / (d * c)
    L[2,0] = -(g2 * v**2 + v * c) / (d * c)
   
    L[0,1] = g1 * v / c**2
    L[1,1] = (c - g1 * v) / (d * c)
    L[2,1] = (c + g1 * v) / (d * c)
   
    L[0,2] = -g1 / c**2
    L[1,2] =  g1 / (d * c)
    L[2,2] = -g1 / (d * c)
   
    # Eigenvector matrix
    R[0,0] = 1.0
    R[1,0] = v
    R[2,0] = 0.5 * v**2
   
    R[0,1] = f
    R[1,1] = (v + c) * f
    R[2,1] = (h + v * c) * f
   
    R[0,2] = -f
    R[1,2] = -(v - c) * f
    R[2,2] = -(h - v * c) * f

    return R, L

# Rusanov flux
def numflux(ul,ur):
    fl0,fl1,fl2 = flux(ul[0],ul[1],ul[2])
    fr0,fr1,fr2 = flux(ur[0],ur[1],ur[2])
    laml = max_eig(ul)
    lamr = max_eig(ur)
    lam  = max(laml,lamr)
    f = np.empty(3)
    f[0] = 0.5*(fl0 + fr0) - 0.5*lam*(ur[0] - ul[0])
    f[1] = 0.5*(fl1 + fr1) - 0.5*lam*(ur[1] - ul[1])
    f[2] = 0.5*(fl2 + fr2) - 0.5*lam*(ur[2] - ul[2])
    return f
