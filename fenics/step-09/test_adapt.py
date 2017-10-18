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
plt.figure()
plot(mesh,interactive=False)
plt.savefig('mesh.png')
# Plot refined initial mesh
plt.figure()
plot(mesh_new,interactive=False)
plt.savefig('mesh_new.png')
plt.show()
