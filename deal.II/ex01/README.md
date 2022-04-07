# Solve BVP in 1-D with Dirichlet bc

Note: A dollar sign indicates the shell prompt you get in a terminal.

## How to run
Generate Makefile
```
$ cmake .
```
Compile the code
```
$ make
```
Run the executable
```
$ ./demo
```
The solution is saved into a file. You can plot it in gnuplot
```
$ gnuplot
gnuplot> load 'plot.gnu'
```
The script `make_eps.gnu` generates solution in eps format, run it as
```
$ gnuplot make_eps.gnu
$ gv sol.eps
```

## Exercise 1

Try to increase degree in the main function
```
int degree = 10;
unsigned int nrefine = 1;
```
In this case, there are just two elements. To get better solution visualization, use more patches, e.g.,
```
data_out.build_patches (2*fe.degree);
```

## Exercise 2
Try a direct solver. Remove the contents of `solve` function and put this
```
SparseDirectUMFPACK solver;
solver.initialize (system_matrix);
solver.vmult (solution, system_rhs);
```
You must include following file
```
#include <deal.II/lac/sparse_direct.h>
```
in order to use the sparse solver.

## Exercise 3
Write a function to the compute error norm in solution and its derivative.

## Exercise 4
Solve the more general problem

-(a u')' + b u' + c u = f

with a(x) = 1 + x, b(x) = x, c(x) = 1. You can keep same exact solution and obtain f(x). Add one class for a, b, c which are derived from `Function<1>` class. Modify the `assemble_system` function. Since the problem is now non-symmetric, we cannot use CG solver. Use the direct solver as explained above.
