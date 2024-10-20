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

## Exercise: output solution in pvtu format

Save solution in vtu format using `write_vtu_with_pvtu_record` function for parallel codes.

## TODO: limiter

Limiters on quad meshes have not been implemented. For some ideas, see

Krishnadutt, Limiting techniques for the discontinuous Galerkin method on unstructured meshes, PhD Thesis.
http://hdl.handle.net/10012/18566

Frederico Bolsoni Oliveira, João Luiz F. Azevedo & Z. J. Wang
ANALYSIS OF THE R FAMILY OF LIMITERS APPLIED TO HIGH-ORDER FR/CPR SCHEMES FOR THE SIMULATION OF SUPERSONIC FLOWS
https://www.icas.org/ICAS_ARCHIVE/ICAS2024/data/papers/ICAS2024_1181_paper.pdf
