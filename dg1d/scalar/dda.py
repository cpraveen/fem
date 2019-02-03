# Plot dissipation/dispersion of physical mode
# Run as
#    python3 ./dda.py
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

def get_eig(degree):
    print('Degree = ', degree)
    k  = degree     # degree
    Nq = k + 1      # number of quadrature points
    nd = k + 1      # number of dofs

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
            # Physical mode is one with max dissipation
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
    return pmode,K,eigr,eigi

for degree in range(1,7):
    pmode,K,eigr,eigi = get_eig(degree)

    # Physical mode: seems to be last
    plt.figure(1)
    plt.plot(K, eigr[:,pmode],lw=2)

    # Physical mode: seems to be last
    plt.figure(2)
    plt.plot(K, eigi[:,pmode],lw=2)


plt.figure(1)
K = np.linspace(0,1,100)
plt.plot(K,K*np.pi,'k--')
plt.ylabel('$\Omega_r/(N+1)$')
plt.xlabel('$K/\pi$')
plt.grid(True)
plt.legend(('N=1','N=2','N=3','N=4','N=5','N=6'))
plt.savefig('omegar_phy.pdf')

plt.figure(2)
plt.ylabel('$\Omega_i/(N+1)$')
plt.xlabel('$K/\pi$')
plt.grid(True)
plt.legend(('N=1','N=2','N=3','N=4','N=5','N=6'))
plt.savefig('omegai_phy.pdf')

plt.show()
