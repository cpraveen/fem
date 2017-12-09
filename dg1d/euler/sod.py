import numpy as np
from gas import *

xmin, xmax = 0.0, 1.0

def initial_condition(x):
    rhol, rhor = 1.0, 0.125
    vell, velr = 0.0, 0.0
    prel, prer = 1.0, 0.1
    rho = np.empty_like(x)
    mom = np.empty_like(x)
    ene = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx < 0.5:
            rho[i] = rhol
            mom[i] = rhol*vell
            ene[i] = prel/(gasGam-1.0) + 0.5*rhol*vell**2
        else:
            rho[i] = rhor
            mom[i] = rhor*velr
            ene[i] = prer/(gasGam-1.0) + 0.5*rhor*velr**2
    return rho,mom,ene
