# H(div) method for Poisson equation

By default, it generates Cartesian grids. You can randomly perturb the grid to get quad elements by setting, e.g.,

```text
set perturb grid = 0.1
```

Increasing this value increases the amount of perturbation.

Observe the convergence rates on Cartesian and quad grids.

> Currently gmres is used without preconditioner since ILU does not work with BlockVector. TODO: Switch to Petsc or Trilinos solvers.

## Dirichlet problem

See `dirichlet.h` file.

```c++
cmake -DPROBLEM=1 .
make release
make
./main
```

## Neumann problem

See `neumann.h` file.

```shell
cmake -DPROBLEM=2 .
make release
make
./main
```
