# DG code for 2-D conservation laws

* `scalar_legendre`: Scalar PDE using Legendre basis on Cartesian grids
* `system_legendre`: System PDE using Legendre basis on Cartesian grids
* `system_legendre_mpi`: System PDE using Legendre basis on Cartesian grids with mpi
* `system_lagrange_mpi`: System PDE using Lagrange basis on Cartesian and quadrilateral (curved) grids with mpi
* `models`: PDE and test cases, use these together with the code in `system_*` directories

## Exercise: Using triangular grids

All the codes are for Cartesian and quadrilateral grids, but deal.II also supports triangular grids now, though it is still under development. Try to write a code using triangles by modifying the `system_legendre_mpi` code. You need to use

* `FE_SimplexDGP`
* `MappingFE`
* `QGaussSimplex`

It looks like barycentric coordinates based on uniformly distributed points are used to define the basis functions. This means that the mass matrix will not be diagonal. You can compute mass matrix on each element as a FullMatrix, invert it and store the inverse mass matrix as a PetscWrappers::SparseMatrix; but then check if you can call PetscWrappers::SparseMatrix.vmult on a deal.II parallel vector. If this is not possible, then you have to use PetscWrappers::Vector for all vectors.

Most other things, except limiting, should remain the same.
