import numpy as np
from gas import *

xmin, xmax = -5.0, 5.0
Tf = 1.8

def initial_condition(x):
    epsilon = 0.2
    velocity = (x < -4.0)*2.629369
    pressure = (x < -4.0)*10.33333 + (x >= -4.0)*1.0
    rho = (x < -4.0)*3.857143 + (x >= -4.0)*(1+epsilon*np.sin(5*x))
    E = pressure/(gasGam-1.0) + 0.5 * rho * velocity**2
    return rho, rho*velocity, E
