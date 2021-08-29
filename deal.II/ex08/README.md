# Laplace equation with discontinuous dirichlet bc: adaptive refinement

This problem was given as an exercise in ex04 where uniform grid refinement is performed.

This code solves Laplace equation

```text
 Laplace(u) = 0  in  (-0.5, +0.5) x (0, 1)
         u  = 1  on  (-0.5,  0.0) x (y = 0)
         u  = 0  on  ( 0.0, +0.5) x (y = 0)
         u  = (1/pi) * atan(y/x)  elsewhere on boundary
```

The exact solution is

```text
u = (1/pi) * atan(y/x)
```

Now compile and run the code

```shell
cmake .
make release
make
./demo
```

## Exercise: Solve on triangular grid

Modify this example to solve on triangular grids. See ex02b for an example of how to do this.

