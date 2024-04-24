# Solve using MatrixCreator

For many standard problems, there are functions available in deal.II to construct the matrix and right hand side, see the namespaces `MatrixCreator` and `VectorTools`.

The code solves

```text
-Laplace(u) = 1 on (0,1)x(0,1)
         u  = 0 on boundary
```

We use `ZeroFunction` for the boundary values and `ConstantFunction` for the right hand side function.

## Exercise: FunctionParser for more general functions

`FunctionParser` is an easy way to define more general functions whose expression is not too complicated. To use this, include

```c++
#include <deal.II/base/function_parser.h>
```

Define right hand side function

```c++
   std::map<std::string,double> constants;
   constants["pi"] = numbers::PI;
   FunctionParser<2> rhs_function(1);
   rhs_function.initialize("x,y",
                           "8*pi^2*sin(2*pi*x)*sin(2*pi*y)",
                           constants);
```

and use `rhs_function` in `VectorTools::create_right_hand_side`. Define the exact solution

```c++
   FunctionParser<2> exact_solution(1);
   exact_solution.initialize("x,y",
                             "sin(2*pi*x)*sin(2*pi*y)",
                             constants);
```

and use `exact_solution` in `VectorTools::interpolate_boundary_values` to apply Dirichlet boundary conditions.

Also, see the documentation

https://www.dealii.org/current/doxygen/deal.II/classFunctionParser.html

## Exercise: Laplace with coefficient function

See the different versions of `MatrixCreator::create_laplace_matrix` function.

## Exercise: Neumann bc

If we have Neumann bc, then the weak formulation has additional terms in the linear functional. These can be assembled using `VectorTools::create_boundary_right_hand_side` function.
