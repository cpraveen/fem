# Laplace equation in non-convex domain: adaptive refinement

This code solves Laplace equation in a Gamma shaped domain on adaptively refined grids. This is similar to `ex06`. The mesh is generated using Gmsh using the geo file in `ex06`

```shell
gmsh -2 ../ex06/Gamma.geo -o Gamma.msh
```

Now compile and run the code

```shell
cmake .
make release
make
./demo
```

This produces a latex file `error.tex` containing error.

```shell
pdflatex error.tex
```

Open and see `error.pdf` for a table of error versus number of dofs. The solution and error are saved in files `solution-00.vtk`, `solution-01.vtk`, etc. Open them together in VisIt and see the change with grid refinement.

```shell
visit -o solution-*.vtk
```

If we compare the error vs dofs from uniform and adaptive refinement, we see that adaptive refinement allows us to achieve same error level as uniform refinement but with substantially fewer dofs.

## Exercise: Solve on triangular grid

Modify this example to solve on triangular grids. See ex02b for an example of how to do this.

## Additional examples

See `ex08`, `ex09`, and `step-6` in the deal.II tutorials

http://www.dealii.org/developer/doxygen/deal.II/step_6.html

for other example on adaptively refined grids.
