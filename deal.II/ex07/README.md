# Laplace equation in non-convex domain: adaptive refinement

This code solves Laplace equation in a Gamma shaped domain on adaptively refined grids. This is similar to `ex06`. The mesh is generated using Gmsh using the geo file in `ex06`
```
gmsh -2 ../ex06/Gamma.geo -o Gamma.msh
```
Now compile and run the code
```
cmake .
make
./demo
```
This produces a latex file `error.tex` containing error and convergence rates.
```
pdflatex error.tex
```
Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate. The solution and error are saved in files `solution-00.vtk`, `solution-01.vtk`, etc. Open them together in VisIt and see the change with grid refinement.
