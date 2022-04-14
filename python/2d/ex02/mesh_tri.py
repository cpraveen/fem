import meshio
import numpy as np
import matplotlib.pyplot as plt

mesh = meshio.read("annulus.msh")

x, y = mesh.points[:,0], mesh.points[:,1]
tri = mesh.cells_dict["triangle"]
faces = mesh.cells_dict["line"]
inner = mesh.cell_sets_dict["inner"]["line"]
outer = mesh.cell_sets_dict["outer"]["line"]

# Show mesh
plt.figure()
plt.triplot(x, y, tri)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Mesh')
plt.axis('equal')
plt.savefig('mesh.svg')

# Show boundary faces
plt.figure()
for f in inner:
    face = faces[f]
    xf = [ x[face[0]], x[face[1]] ]
    yf = [ y[face[0]], y[face[1]] ]
    plt.plot(xf, yf, '-r')
for f in outer:
    face = faces[f]
    xf = [ x[face[0]], x[face[1]] ]
    yf = [ y[face[0]], y[face[1]] ]
    plt.plot(xf, yf, '-b')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Boundary faces')
plt.axis('equal')

plt.show()
