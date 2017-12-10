from riemann import *

xmin, xmax = 0.0, 1.0
Tf = 0.2

def initial_condition(x):
    rhol, rhor = 1.0, 0.125
    vell, velr = 0.0, 0.0
    prel, prer = 1.0, 0.1
    xs = 0.5
    vl = [rhol, vell, prel]
    vr = [rhor, velr, prer]
    return riemann_data(vl,vr,xs,x)
