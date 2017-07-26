"""
Test grid adaptation
"""
from dolfin import *
import matplotlib.pyplot as plt

# Initial mesh
n = 10
mesh = UnitSquareMesh(n,n)

flag = CellFunction("bool", mesh)
for c in cells(mesh):
   if c.index() == 50:
      flag[c] = True
   else:
      flag[c] = False

# Refine the mesh
mesh_new = refine(mesh, flag)

# Plot initial mesh
plot(mesh,interactive=False)
plt.show()
plt.savefig('mesh.png')
# Plot refined initial mesh
plot(mesh_new,interactive=False)
plt.show()
plt.savefig('mesh_new.png')
