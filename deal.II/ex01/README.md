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

Write a function to the compute error norm in solution and its derivative. The following function does this for the solution error. We use the `BoundaryValues` class to compute the exact solution.

```c++
void LaplaceProblem::compute_error () const
{
   std::cout << "Computing error norm ...\n";
   QGauss<1>  quadrature_formula(fe.degree + 3);
   FEValues<1> fe_values (fe, quadrature_formula,
                          update_values   | update_gradients |
                          update_quadrature_points | update_JxW_values);

   const unsigned int   n_q_points    = quadrature_formula.size();

   std::vector<double> ex_sol_values (n_q_points);
   std::vector<double> num_sol_values (n_q_points);
   const BoundaryValues exact_solution;

   double u_err = 0.0;
   for(const auto &cell : dof_handler.active_cell_iterators())
   {
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, num_sol_values);
      exact_solution.value_list(fe_values.get_quadrature_points(),
                                ex_sol_values);
      double cell_u_err = 0.0;
      for(const auto q : fe_values.quadrature_point_indices())
      {
         double err = ex_sol_values[q] - num_sol_values[q];
         cell_u_err += pow(err, 2) * fe_values.JxW(q);
      }
      u_err += cell_u_err;
   }
   u_err = sqrt(u_err);
   std::cout << "L2 error norm = " << u_err << std::endl;
}
```

Run with `nrefine = 5` and then with `nrefine = 6`, and use the two error norm values to estimate convergence rate.

## Exercise 4

Solve the more general problem

```
-(a u')' + b u' + c u = f
```

with `a(x) = 1 + x, b(x) = x, c(x) = 1`. You can keep same exact solution and obtain `f(x)`. Add one class for a, b, c which are derived from `Function<1>` class. Modify the `assemble_system` function. Since the problem is now non-symmetric, we cannot use CG solver. Use the direct solver as explained above.
