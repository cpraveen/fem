# Linear advection in 2d

We assume velocity is divergence-free so we can write in conservation form.

## rotate.h: rotating gaussian profile

Solve in square domain, can be used with legendre or lagrange DG code.

```text
set mapping = cartesian
set grid    = 100,100
```

You can use the code in `system_legendre` or `system_legendre_mpi` or `system_lagrange_mpi` for this problem.

## rotate_annulus.h: rotating gaussian profile

Solve in annular domain, can be used with lagrange DG code.

Generate mesh

```shell
gmsh -2 annulus.geo
```

and set parameters

```text
set basis   = gl
set mapping = q
set grid    = annulus.msh
```

You can use the code in `system_lagrange_mpi` for this problem.
