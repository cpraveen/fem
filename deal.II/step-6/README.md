# 2-D Poisson equation with discontinuous coefficient

This is based on step-6 from the deal.II tutorials but modified to use grid made in gmsh. The changes are in

* run: Reads gmsh grid and attaches SphericalManifold
* output_results: enable curved vtk output (use Paraview)

```shell
cmake .
make release && make
gmsh -2 ../grid/circ_disc_nice.geo -o circ_disc_nice.msh
./step-6
```
