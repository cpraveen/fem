# 2d system of conservation law

Solves 2d system of conservation laws on Cartesian grids using Legendre basis. The domain need not be a rectangle, but can be anything, but the cells need to be rectangles.

This code is very similar to the one in `system_legendre` but this one is parallel. See the differences

```shell
code -n -d main.cc ../system_legendre/main.cc
code -n -d dg.h    ../system_legendre/dg.h
```

## Run an example

```shell
ln -s ../models/euler/pde.h
ln -s ../models/euler/isentropic_vortex.h problem.h
cmake .
make release
make
mpirun -np 4 ./main > log.txt 2>&1 &
tail -f log.txt
visit -o solution.xdmf
```
