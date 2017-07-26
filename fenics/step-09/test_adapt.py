"""
Test grid adaptation
"""
from dolfin import *

# Initial mesh
n = 10
mesh = UnitSquareMesh(n,n)

flag = CellFunction("bool", mesh)
for c in cells(mesh):
   if c.index() == 50:
      flag[c] = True
   else:
      flag[c] = False

mesh_new = refine(mesh, flag)
File('mesh.pvd') << mesh
File('mesh_new.pvd') << mesh_new
