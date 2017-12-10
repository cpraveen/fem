import numpy as np
from gas import *

def riemann_data(vl,vr,xs,x):
    rhol, rhor = vl[0], vr[0]
    vell, velr = vl[1], vr[1]
    prel, prer = vl[2], vr[2]
    rho = np.empty_like(x)
    mom = np.empty_like(x)
    ene = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx < xs:
            rho[i] = rhol
            mom[i] = rhol*vell
            ene[i] = prel/(gasGam-1.0) + 0.5*rhol*vell**2
        else:
            rho[i] = rhor
            mom[i] = rhor*velr
            ene[i] = prer/(gasGam-1.0) + 0.5*rhor*velr**2
    return rho,mom,ene
