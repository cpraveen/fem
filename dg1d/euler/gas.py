import numpy as np

gasGam = 1.4
gasR   = 1.0

def pressure(rho,mom,ene):
    return (gasGam-1.0)*(ene - 0.5*mom**2/rho)

def flux(rho,mom,ene):
    pre  = pressure(rho,mom,ene)
    frho = mom
    fmom = pre + mom**2/rho
    fene = (ene + pre)*mom/rho
    return frho,fmom,fene

def max_eig(u):
    vel = u[1]/u[0]
    pre = pressure(u[0],u[1],u[2])
    lam = np.abs(vel) + np.sqrt(gasGam * pre / u[0])
    return lam

# Based on cell average value
def max_speed(rho,mom,ene):
    vel = mom[:,0] / rho[:,0]
    pre = pressure(rho[:,0], mom[:,0], ene[:,0])
    lam = np.abs(vel) + np.sqrt(gasGam * pre / rho[:,0])
    return lam.max()

def numflux(ul,ur):
    fl0,fl1,fl2 = flux(ul[0],ul[1],ul[2])
    fr0,fr1,fr2 = flux(ur[0],ur[1],ur[2])
    laml = max_eig(ul)
    lamr = max_eig(ur)
    lam  = max(laml,lamr)
    f = np.zeros(3)
    f[0] = 0.5*(fl0 + fr0) - 0.5*lam*(ur[0] - ul[0])
    f[1] = 0.5*(fl1 + fr1) - 0.5*lam*(ur[1] - ul[1])
    f[2] = 0.5*(fl2 + fr2) - 0.5*lam*(ur[2] - ul[2])
    return f
