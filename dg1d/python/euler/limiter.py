import numpy as np

def minmod(a,b,c,Mdx2):
    if np.abs(a) < Mdx2:
        return a
    sa = np.sign(a)
    sb = np.sign(b)
    sc = np.sign(c)
    if sa==sb and sb==sc:
        return sa * np.abs([a,b,c]).min()
    else:
        return 0.0
