import numpy as np
from gas import *

def riemann_data(vl,vr,xs,x):
    rhol, rhor = vl[0], vr[0]
    vell, velr = vl[1], vr[1]
    prel, prer = vl[2], vr[2]
    enel = prel/(gasGam-1.0) + 0.5*rhol*vell**2
    ener = prer/(gasGam-1.0) + 0.5*rhor*velr**2
    rho = (x <= xs)*rhol      + (x > xs)*rhor
    mom = (x <= xs)*rhol*vell + (x > xs)*rhor*velr
    ene = (x <= xs)*enel      + (x > xs)*ener
    return rho,mom,ene
