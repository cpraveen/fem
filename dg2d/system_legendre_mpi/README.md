# 2d system of conservation law on Cartesian grids

Solves 2d system of conservation laws on Cartesian grids using Legendre basis. The domain need not be a rectangle, but can be anything, but the cells need to be rectangles whose sides are aligned with the x, y axes.

This code is very similar to the one in `system_legendre` but this one is parallel. See the differences

```shell
code -n -d main.cc ../system_legendre/main.cc
code -n -d dg.h    ../system_legendre/dg.h
```

## Run an example

In a terminal, go to `fem/dg2d/system_legendre_mpi`

```shell
ln -s ../models/euler/pde.h
ln -s ../models/euler/isentropic_vortex/problem.h
cmake .
make release
make
```

In another terminal, go to `fem/dg2d/models/euler/isentropic_vortex`

```shell
mpirun -np 4 ../../../system_legendre_mpi/main input.prm > log.txt 2>&1 &
tail -f log.txt
visit -o solution.xdmf
```

When the code is running, if you use `top`, you should see four instances of `main` program running.
