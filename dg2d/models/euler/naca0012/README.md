# Flow over NACA0012 airfoil

Set angle of attack and mach number in `problem.h` file.

```shell
gmsh -2 naca.geo
mpirun -np 4 ../../../system_lagrange_mpi/main input.prm > log.txt 2>&1 &
visit -o solution.xdmf
```
