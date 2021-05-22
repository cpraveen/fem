# Laplace equation in non-convex domain

This code solves Laplace equation in a Gamma shaped domain on succesively refined grids. The mesh is generated using Gmsh
```
gmsh -2 Gamma.geo
```
Now compile and run the code
```
cmake .
make release
make
./demo
```
This produces a latex file `error.tex` containing error and convergence rates.
```
pdflatex error.tex
```
Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate. The solution and error are saved in files `solution-00.vtk`, `solution-01.vtk`, etc. Open them together in VisIt and see the change with grid refinement.
