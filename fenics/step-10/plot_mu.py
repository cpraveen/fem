'''
Visualize the discontinuous coefficient in the Poisson equation

Before running this code, generate mesh like this:
    gmsh -2 -format msh2 mesh.geo
    dolfin-convert -i gmsh -o xml mesh.msh mesh.xml
'''
from dolfin import *
import matplotlib.pyplot as plt

mu1 = 1.0
mu2 = 10.0
mu = Expression('x[0]*x[0] + x[1]*x[1] <= 0.5*0.5 ? mu1 : mu2',
                mu1=mu1, mu2=mu2, degree=0)

def plot_mu(mesh,fname):
    V0 = FunctionSpace(mesh,'DG',0)
    mu0 = project(mu,V0)

    plt.figure(figsize=(12,6))
    plt.subplot(121)
    plot(mesh)
    plt.subplot(122)
    c = plot(mu0)
    plt.colorbar(c)
    plt.savefig(fname)
    print('Coefficient mu printed to '+fname)

n = 10
mesh = RectangleMesh(Point(-1.0, -1.0), Point(+1.0, +1.0), n, n)
plot_mu(mesh,'mu1.pdf')

mesh = Mesh('mesh.xml')
plot_mu(mesh,'mu2.pdf')
