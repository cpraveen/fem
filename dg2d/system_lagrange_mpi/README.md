# 2d system of conservation law on unstructured grids

Solves 2d system of conservation laws on unstructured grids using Lagrange basis.

## Run an example

Set parameters in `input.prm` file.

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

## Exercise: Flow over cylinder (euler)

Solve subsonic flow over cylinder at Mach number of 0.3; make a grid in Gmsh and run the code for a long time to reach steady solution.

## Exercise: Ringleb flow (euler)

The problem is described in Section 7.13.6 of this book

http://www.ossanworld.com/cfdbooks/download_idolikecfd_free_2nd_edition.html

Generate a mesh for this problem using Gmsh and set up the problem file. Run it for a long enough time that steady solution is reached.
