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

## Compile and run

```shell
cmake .
make release
make
mpirun -np 4 ../../../system_lagrange_mpi/main input.prm > log.txt 2>&1 &
```

## Viewing results

Using visit

```shell
visit -o ./solution.xdmf
```

Or use visit python script

```shell
visit -cli -s ./plot_visit.py
```

Using pyvista

```shell
python plot_pyvista.py
```
