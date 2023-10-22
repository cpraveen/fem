# Fourier stability analysis of RKDG scheme
import numpy as np
from numpy.linalg import eigvals
import matplotlib.pyplot as plt
import argparse
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

nwave = nd*500
wavenums = np.linspace(0,nd*np.pi,nwave)
eigr = np.zeros((nwave,nd))
eigi = np.zeros((nwave,nd))

for i,kdx in enumerate(wavenums):
    C = A + np.exp(-1j*kdx)*Bm - Bp
    eig = 1j * eigvals(C)
    if i == 0:
        eigr[i,:] = np.real(eig)
        eigi[i,:] = np.imag(eig)
        # Physical mode has zero dissipation
        pmode = np.argmax(eigi[i,:])
    else:
        # Find closest eigenvalue to previous one
        eig_old = eigr[i-1,:] + 1j * eigi[i-1,:]
        for j in range(nd):
            jj = np.argmin( np.abs(eig_old[j] - eig) )
            eigr[i,j] = np.real(eig[jj])
            eigi[i,j] = np.imag(eig[jj])

K = wavenums/nd/np.pi
eigr = eigr/nd
eigi = eigi/nd

# Physical mode
plt.figure()
plt.plot(K, eigr[:,pmode],lw=2)
plt.plot(K, K*np.pi, 'k--')
plt.ylabel('$\Omega_r/(N+1)$')
plt.xlabel('$K/\pi$')
plt.title('Dispersion: Physical mode, Degree, N = '+str(k))
plt.grid(True)

# Physical mode
plt.figure()
plt.plot(K, eigi[:,pmode],lw=2)
plt.ylabel('$\Omega_i/(N+1)$')
plt.xlabel('$K/\pi$')
plt.title('Dissipation: Physical mode, Degree, N = '+str(k))
plt.grid(True)

# all modes: real
plt.figure()
for i in range(nd):
    plt.plot(K, eigr[:,i],lw=2)
plt.plot(K, K*np.pi, 'k--')
plt.ylabel('$\Omega_r/(N+1)$')
plt.xlabel('$K/\pi$')
plt.title('Dispersion: All modes, Degree, N = '+str(k))
plt.grid(True)
plt.savefig('omegar_all.pdf')

# all modes: imag
plt.figure()
for i in range(nd):
    plt.plot(K, eigi[:,i],lw=2)
plt.ylabel('$\Omega_i/(N+1)$')
plt.xlabel('$K/\pi$')
plt.title('Dissipation: All modes, Degree, N = '+str(k))
plt.grid(True)
plt.savefig('omegai_all.pdf')

plt.show()
