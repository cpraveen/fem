# Solve Poisson equation in 2d/3d

Steps to run this example

```
cmake .
make
./demo
```

## Run the problem in 3d

To run in 3-d, change the problem in `main` function as

```
   LaplaceProblem<3> problem (degree);
```

save, recompile and run the code.

## Exercise: Zero Neumann bc

Let us solve the problem with same right hand side, but apply Dirichlet bc only on left and right side of the domain. Modify the grid generation as

```
   GridGenerator::hyper_cube (triangulation, 0, 1, true);
```

This assigns different boundary indicators for the faces on the four sides of the square (0=left, 1=right, 2=bottom, 3=top). To apply Dirichlet bc on left/right, call `VectorTools::interpolate_boundary_values` once with boundary indicator = 0 and once with boundary indicator = 1 and then call `MatrixTools::apply_boundary_values`. Nothing else needs to be done to apply zero Neumann bc. Note that solution for this case will be different compared to the pure Dirichlet case.
