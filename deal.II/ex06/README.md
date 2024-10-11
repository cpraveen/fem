# Laplace equation in non-convex domain

This code solves Laplace equation in a Gamma shaped domain on succesively refined grids. The mesh is generated using Gmsh

```shell
gmsh -2 Gamma.geo
```

Now compile and run the code

```shell
cmake .
make release
make
./demo
```

This produces a latex file `error.tex` containing error and convergence rates.

```shell
pdflatex error.tex
```

Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate. The solution and error are saved in files `solution-00.vtu`, `solution-01.vtu`, etc. Open them together in VisIt and see the change with grid refinement.

```shell
visit -o solution-*.vtu
```

## Exercise: Solve on triangular grid

Modify this example to solve on triangular grids. See ex02b for an example of how to do this.
