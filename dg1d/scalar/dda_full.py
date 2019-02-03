# Fourier stability analysis of RKDG scheme
import numpy as np
from numpy.linalg import eigvals
import argparse
import matplotlib.pyplot as plt
from basis import *
from matplotlib import rcParams
rcParams['font.size'] = 12
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 2
rcParams['lines.markersize'] = 2
rcParams['axes.titlesize'] = 12
rcParams['axes.labelsize'] = 12

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-degree', type=int, help='Degree', required=True)
parser.add_argument('-cfl', type=float, help='CFL', required=True)
parser.add_argument('-scheme',
                choices=('fe','ssprk22','ssprk33','ssprk43','ssprk54','rk4','ts2'),
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

def amplification_matrix_two_stage(scheme, nu, C1, C2):
    if scheme == 'ts2':
        a21 = 0.50
        A21 = 1.0/8.0

        b1 = 1.0
        b2 = 0.0
        B1 = 1.0/6.0
        B2 = 1.0/3.0

        G1 = I + a21*nu*C1 - 2.0*A21*nu*nu*C2
        G2 = I + b1*nu*C1 - 2.0*B1*nu*nu*C2
        G3 = b2*nu*C1 - 2.0*B2*nu*nu*C2
        H = G2 + G3@G1
    else:
        print('Unknown time scheme')
        exit()
    return H

M = np.zeros((nd,nd))
A_1 = np.zeros((nd,nd))
A_2 = np.zeros((nd,nd))
B1m= np.zeros((nd,nd))
B1p= np.zeros((nd,nd))
B2m= np.zeros((nd,nd))
B2p= np.zeros((nd,nd))
for i in range(nd):
    for j in range(nd):
        B1m[i,j] = shape_value(i,-1.0) * shape_value(j,+1.0)
        B1p[i,j] = shape_value(i,+1.0) * shape_value(j,+1.0)

        B2m[i,j] = shape_value(i,-1.0) * shape_grad(j,+1.0)  
        B2p[i,j] = shape_value(i,+1.0) * shape_grad(j,+1.0)  
        for q in range(Nq):
            M[i,j] += 0.5*Vf[q,i]*Vf[q,j]*wg[q]
            A_1[i,j] += Vg[q,i]*Vf[q,j]*wg[q]
            A_2[i,j] += Vg[q,i]*Vg[q,j]*wg[q]
print("M=",M)
print("A1=",A_1)
print("A2=",A_2)

nu = args.cfl
nwave    = nd*500
wavenums = np.linspace(0,nd*np.pi,nwave)
eigr = np.zeros((nwave,nd))
eigi = np.zeros((nwave,nd))

for i,kdx in enumerate(wavenums):
    if args.scheme == 'ts2':
        C1 = A_1 + np.exp(-1j*kdx)*B1m - B1p
        C2 = A_2 + np.exp(-1j*kdx)*B2m - B2p
        H = amplification_matrix_two_stage(args.scheme, nu, C1, C2)
    else:
        C= A_1 + np.exp(-1j*kdx)*B1m - B1p
        H = amplification_matrix(args.scheme, nu, C)
    eig = eigvals(H)
    if i == 0:
        eigr[i,:] = np.real(eig)
        eigi[i,:] = np.imag(eig)
        # Physical has least dissipation at zero wavenumber
        pmode = np.argmax( np.abs(eig) )
    else:
        # Find closest eigenvalue to previous one
        eig_old = eigr[i-1,:] + 1j * eigi[i-1,:]
        for j in range(nd):
            jj = np.argmin( np.abs(eig_old[j] - eig) )
            eigr[i,j] = np.real(eig[jj])
            eigi[i,j] = np.imag(eig[jj])

# String for plot title
tstr = ', N = '+str(k)+', scheme = '+args.scheme

eigv = eigr + 1j * eigi
diss =  np.abs  (eigv)
disp = -np.angle(eigv)/np.pi
K    = wavenums/np.pi/nd

# physical mode: usually last one, but check
plt.figure()
plt.plot(K, diss[:,pmode], lw=2)
plt.ylabel('Magnitude of eigenvalue')
plt.xlabel('$K/\pi$')
plt.title('Dissipation: Physical mode'+tstr)
plt.grid(True)
plt.savefig('diss_phy_deg'+str(k)+'_'+args.scheme+'_cfl'+str(round(nu,3))+'.pdf')

plt.figure()
plt.plot(K, disp[:,pmode],lw=2)
plt.plot(K, wavenums*nu/np.pi, 'k--', lw=2)
plt.ylabel('Angle/$\pi$ of eigenvalue')
plt.xlabel('$K/\pi$')
plt.title('Dissipation: Physical mode'+tstr)
plt.grid(True)
plt.savefig('disp_phy_deg'+str(k)+'_'+args.scheme+'_cfl'+str(round(nu,3))+'.pdf')

# all modes: dissipation
plt.figure()
for i in range(nd):
    plt.plot(K, diss[:,i], lw=2)
plt.ylabel('Magnitude of eigenvalue')
plt.xlabel('$K/\pi$')
plt.title('Dissipation: all modes'+tstr)
plt.grid(True)
plt.savefig('diss_all_deg'+str(k)+'_'+args.scheme+'_cfl'+str(round(nu,3))+'.pdf')

# all modes: dispersion
plt.figure()
for i in range(nd):
    plt.plot(K, disp[:,i], lw=2)
plt.plot(K, wavenums*nu/np.pi, 'k--', lw=2)
plt.ylabel('Angle/$\pi$ of eigenvalue')
plt.xlabel('$K/\pi$')
plt.title('Dispersion: all modes'+tstr)
plt.grid(True)
plt.savefig('disp_all_deg'+str(k)+'_'+args.scheme+'_cfl'+str(round(nu,3))+'.pdf')

plt.show()
