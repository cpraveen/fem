# H(div) method for Poisson equation

Laplace equation

$$
-\Delta p = f
$$

is written as first order system

$$
J = \nabla p, \qquad \nabla \cdot J = -f
$$

**Dirichlet problem**: weakly imposed

$$
p = 0 \quad \textrm{on} \quad \partial\Omega
$$

**Neumann problem**: strongly imposed

$$
J \cdot n = 0 \quad \textrm{on} \quad \partial\Omega
$$

together with constraint

$$
\int_{\partial\Omega}p ds = 0
$$

**Weak formulation**: in either case it is

$$
(J, \tilde{J}) + (\nabla\cdot\tilde{J},p) = 0, \qquad
(\nabla\cdot J, \tilde{p}) = -(f, \tilde{p})
$$

The FE spaces are: RT(k) for $J$ and DG(k) for $p$ where $k \ge 0$. On Cartesian grids, RT(k) = Q(k+1,k) x Q(k,k+1).

By default, the code generates Cartesian grids. You can randomly perturb the grid to get quad elements by setting, e.g.,

```text
set perturb grid = 0.1
```

Increasing this value increases the amount of perturbation.

Observe the convergence rates on Cartesian and quad grids. Try degree = 0 and degree = 1 on both Cartesian and quad grids. Note that with degree = 0, the divergence error does not converge on quad meshes. This was investigated by Bochev and Ridzal.

> The Schur solver is based on [step-20](https://dealii.org/current/doxygen/deal.II/step_20.html), see step-20 documentation for a discussion of its limitations. Currently gmres is used without preconditioner since ILU does not work with BlockVector. TODO: Switch to Petsc or Trilinos solvers.

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

## References

1. Bochev and Ridzal, Rehabilitation of the lowest-order Raviart–Thomas element on quadrilateral grids, https://doi.org/10.1137/070704265
