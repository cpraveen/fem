# Advection of isentropic vortex

The vortex moves with constant velocity and we use periodic bc.

Set parameters in `input.prm` file.

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

This geo file can also generate a cartesian grid, see comments inside.

## Compile and run

In one terminal

```shell
cd fem/dg2d/system_lagrange_mpi
ln -s ../models/euler/pde.h
ln -s ../models/euler/isentropic_vortex/problem.h
cmake .
make release
make
mpirun -np 4 ../../../system_lagrange_mpi/main input.prm > log.txt 2>&1 &
```

In another terminal

```shell
cd fem/dg2d/models/euler/isentropic_vortex
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

Use visit to create png images at all times

```shell
rm -f density*.png
visit -cli -nowin -s ./make_png.py
open density*.png
```

Plot using pyvista

```shell
python plot_pyvista.py
```
