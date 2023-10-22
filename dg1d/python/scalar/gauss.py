import numpy as np

xmin, xmax = 0.0, 1.0

def initial_condition(x):
    return np.exp(-100*(x-0.5)**2)


