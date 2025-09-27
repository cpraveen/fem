# Parallel fully distributed triangulation

Shows how to create a p::f::t triangulation.

```shell
cmake .
make release
make
mpirun -np 4 ./demo
visit -o grid_000.pvtu
```

Also see the `log_#.txt` files.
