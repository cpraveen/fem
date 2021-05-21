# Poisson with discontinuous coefficient

Run the code by setting ```refine_type = 'uniform'``` and  ```refine_type = 'adaptive'``` and visualize the meshes.

On a general mesh, the discontinuous coefficient cannot be represented properly. There is a mesh created in Gmsh which can resolve the coefficient better. We can visualize it as follows
```
gmsh -2 -format msh2 mesh.geo
dolfin-convert -i gmsh -o xml mesh.msh mesh.xml
python3 plot_mu.pdf
open mu1.pdf mu2.pdf
```

## Exercise: Use the Gmsh mesh and run the adaptive refinement
Read mesh like this
```
mesh = Mesh('mesh.xml')
```
