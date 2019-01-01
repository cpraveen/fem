# Fourier stability analysis of RKDG scheme
import numpy as np
from numpy.linalg import eigvals
import argparse
from basis import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-degree', type=int, help='Degree', required=True)
parser.add_argument('-cfl_min', type=float, help='Min cfl', default=0.0)
parser.add_argument('-cfl_max', type=float, help='Max cfl', default=1.0)
parser.add_argument('-scheme', choices=('ssprk22','ssprk33'),
                    help='Time scheme', required=True)
args = parser.parse_args()

k  = args.degree     # degree
Nq = k + 1           # number of quadrature points
nd = k + 1           # number of dofs

# QGauss position and weights
xg, wg = np.polynomial.legendre.leggauss(Nq)

# Construct Vandermonde matrix for gauss points
Vf = np.zeros((Nq,nd))
Vg = np.zeros((Nq,nd))
for i in range(Nq):
    for j in range(nd):
        Vf[i,j] = shape_value(j, xg[i])
        Vg[i,j] = shape_grad (j, xg[i])

# Identity
I   = np.eye(nd)

def amplification_matrix(scheme, G):
    if scheme == 'ssprk22':
        H = 0.5*(I + G@G)
    elif scheme == 'ssprk33':
        H = (1.0/3.0)*I + (1.0/2.0)*G + (1.0/6.0)*(G@G@G)
    else:
        print('Unknown time scheme')
        exit()
    return H

M = np.zeros((nd,nd))
A = np.zeros((nd,nd))
Bm= np.zeros((nd,nd))
Bp= np.zeros((nd,nd))
for i in range(nd):
    for j in range(nd):
        Bm[i,j] = shape_value(i,-1.0) * shape_value(j,+1.0)
        Bp[i,j] = shape_value(i,+1.0) * shape_value(j,+1.0)
        for q in range(Nq):
            M[i,j] += 0.5*Vf[q,i]*Vf[q,j]*wg[q]
            A[i,j] += Vg[q,i]*Vf[q,j]*wg[q]

print("M=",M)
print("A=",A)

wavenums = np.linspace(0,2*np.pi,500)
cfls  = np.linspace(args.cfl_min,args.cfl_max,100)
for nu in cfls:
    maxeig = 0.0
    for kdx in wavenums:
        C = A + np.exp(-1j*kdx)*Bm - Bp
        G = I + nu*C
        H = amplification_matrix(args.scheme, G)
        eig = np.abs(eigvals(H)).max()
        if eig > maxeig:
            maxeig = eig
    print(nu,maxeig)
    if maxeig > 1.0:
        break
