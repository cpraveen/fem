from riemann import *

xmin, xmax = -5.0, 5.0
Tf = 1.3

def initial_condition(x):
    rhol, rhor = 0.445, 0.5
    vell, velr = 0.698, 0.0
    prel, prer = 3.528, 0.571
    xs = 0.0
    vl = [rhol, vell, prel]
    vr = [rhor, velr, prer]
    return riemann_data(vl,vr,xs,x)
