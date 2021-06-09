import meshio
import numpy as np
import matplotlib.pyplot as plt

mesh = meshio.read("mesh_tri.msh")

x, y = mesh.points[:,0], mesh.points[:,1]
tri = mesh.cells_dict["triangle"]
faces = mesh.cells_dict["line"]

# Show mesh
plt.figure()
plt.triplot(x, y, tri)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Mesh')
plt.axis('equal')

# Show boundary faces
plt.figure()
for face in faces:
    xf = [ x[face[0]], x[face[1]] ]
    yf = [ y[face[0]], y[face[1]] ]
    plt.plot(xf, yf)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Boundary faces')
plt.axis('equal')

# Some smooth function of (x,y)
u = (0.5/(2.0*np.pi)**2) * np.sin(2.0*np.pi*x) * np.cos(2.0*np.pi*y)

plt.figure()
c = plt.tricontourf(x, y, tri, u, cmap='rainbow', levels=20)
plt.colorbar(c)
plt.tricontour(x, y, tri, u, colors='k', levels=20)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Solution')
plt.axis('equal')

plt.show()
