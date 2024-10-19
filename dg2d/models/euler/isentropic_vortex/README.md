# Advection of isentropic vortex

The vortex moves with constant velocity and we use periodic bc.

To run with Cartesian grid, set

```text
set mapping = cartesian
set grid    = 100,100
```

To run with unstructured grid, set

```text
set mapping = q,1
set grid    = grid.msh
```

and generate grid

```shell
gmsh -2 grid.geo
```
