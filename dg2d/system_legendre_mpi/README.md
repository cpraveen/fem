# 2d system of conservation law

Solves 2d system of conservation laws on Cartesian grids using Legendre basis. The domain need not be a rectangle, but can be anything, but the cells need to be rectangles.

## Run an example

```shell
ln -s ../system_legendre/euler/pde.h
ln -s ../system_legendre/euler/isentropic_vortex.h problem.h
cmake .
make release
make
mpirun -np 4 ./main > log.txt 2>&1 &
tail -f log.txt
visit -o solution.xdmf
```
