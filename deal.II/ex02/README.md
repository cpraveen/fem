# Solve Poisson equation in 2d/3d

To run in 3-d, change the problem in main function as
```
   LaplaceProblem<3> problem (degree);
```
save, recompile and run the code.

## Zero Neumann bc 
Let us solve the problem with same right hand side, but apply Dirichlet bc only on left and right side of the domain. Modify the grid generation as
```
   GridGenerator::hyper_cube (triangulation, 0, 1, true);
```
This assigns different boundary indicators for the faces on the four sides of the square (0=left, 1=right, 2=bottom, 3=top). To apply Dirichlet bc on left/right, add an extra call to `VectorTools::interpolate_boundary_values` and `MatrixTools::apply_boundary_values` with boundary indicator 1. Nothing else needs to be done to apply zero Neumann bc.
