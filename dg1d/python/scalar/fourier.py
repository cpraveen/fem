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
parser.add_argument('-scheme',
                choices=('fe','ssprk22','ssprk33','ssprk43','ssprk54','rk4'),
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

def amplification_matrix(scheme, nu, C):
    if scheme == 'fe':
        H = I + nu*C
    elif scheme == 'ssprk22':
        G = I + nu*C
        H = 0.5*(I + G@G)
    elif scheme == 'ssprk33':
        G = I + nu*C
        H = (1.0/3.0)*I + (1.0/2.0)*G + (1.0/6.0)*(G@G@G)
    elif scheme == 'ssprk43':
        G = I + 0.5*nu*C
        H = (2.0/3.0)*G + (1.0/3.0)*G@G@G@G
    elif scheme == 'ssprk54':
        c11 = 0.391752226571890

        c21 = 0.444370493651235
        c22 = 1.0 - c21
        c23 = 0.368410593050371

        c31 = 0.620101851488403
        c32 = 1.0 - c31
        c33 = 0.251891774271694

        c41 = 0.178079954393132
        c42 = 1.0 - c41
        c43 = 0.544974750228521

        c51 = 0.517231671970585
        c52 = 0.096059710526147
        c53 = 1.0 - (c51 + c52)
        c54 = 0.063692468666290
        c55 = 0.226007483236906

        G1 = I + c11*nu*C
        G2 = c21*I + (c22*I + c23*nu*C) @ G1
        G3 = c31*I + (c32*I + c33*nu*C) @ G2
        G4 = c41*I + (c42*I + c43*nu*C) @ G3
        G5 = c52*I + c54*nu*C
        G6 = c53*I + c55*nu*C

        H = c51*G2 + G5@G3 + G6@G4
    elif scheme == 'rk4':
        G1 = nu*C
        G2 = G1 + 0.5*G1@G1
        G3 = G1 + 0.5*G1@G2
        G4 = G1 +     G1@G3
        H =   I + (1.0/6.0)*(G1+G4) + (1.0/3.0)*(G2+G3)
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
        H = amplification_matrix(args.scheme, nu, C)
        eig = np.abs(eigvals(H)).max()
        if eig > maxeig:
            maxeig = eig
    print(nu,maxeig)
    if maxeig - 1.0 > 1.0e-12:
        break
