# Error convergence under uniform refinement

This code solves Poisson equation on unit square on succesively refined grids.
```
cmake .
make
./demo
```
This produces a latex file `error.tex`
```
pdflatex error.tex
```
Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate.

The solution error as a function of x,y is saved into different files as the grid is refined. You can open the error-##.vtk files in visit and see the error distribution in the domain.
