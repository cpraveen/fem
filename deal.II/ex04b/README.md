# Solve Poisson equation in parallel

This is same test case as ex04. We solve it in parallel using `parallel::distributed::Triangulation` and PETSc. Diff the source files to see the differences.

```shell
cmake .
make release
make
mpirun -np 2 ./demo
visit -o sol.pvd
```

You can see the mesh partitions in Visit by `Add -> Mesh -> mesh` and `Add -> Subset -> blocks` and then click `Draw`.
