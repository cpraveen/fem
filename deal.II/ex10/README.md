# Poisson equation with curved boundaries

This code solves Poisson equation on unit square on succesively refined grids.

```shell
gmsh -2 annulus.geo
cmake .
make release && make
./demo
```

This produces a latex file `error.tex`

```shell
pdflatex error.tex
```

Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate.

Run with different degree for solution and mapping.

```
degree = 1; mapping_degree = 1;
degree = 2; mapping_degree = 1;
degree = 2; mapping_degree = 2;
```
