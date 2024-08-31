# 1-d linear scalar PDE using DG

This uses Legendre polynomials as basis functions.

To solve constant linear advection equation

```shell
cp linadv/*.h .
```

To solve Burgers equation

```shell
cp burgers/*.h .
```

Create Makefile (this only needs to be done once)

```shell
cmake .
make release
```

Compile

```shell
make
```

Set parameters in file `input.prm` and run

```shell
rm -f *.gpl
./dg input.prm
gnuplot anim.gnu
```

If you modify `dg.cc`, then you must compile and run again.

## Grid convergence test

Set `nrefine > 1` in `input.prm` and run the code to perform error convergence study. A tex file `error.tex` with convergence table will also be generated. You can compile and see the pdf.

```shell
pdflatex error.tex
open error.pdf
```

## Other ways to specify test data

We have derived classes from `Function` to define initial condition and exact solution.

### Using FunctionParser

Add include file

```c++
#include <deal.II/base/function_parser.h>
```

Define functions

```c++
FunctionParser<1> initial_condition("sin(2*pi*x)");
FunctionParser<1> exact_solution("sin(2*pi*(x-t)");
exact_solution.set_time(param.final_time);
param.xmin = 0.0; param.xmax = 1.0;
```

The gradient is computed using finite differences, so do not depend on it.

### Using SymbolicFunction

Add include file

```c++
#include <deal.II/base/symbolic_function.h>
```

Define functions

```c++
Functions::SymbolicFunction<1> initial_condition("sin(2*pi*x)");
Functions::SymbolicFunction<1> exact_solution("sin(2*pi*(x-t)");
exact_solution.set_time(param.final_time);
param.xmin = 0.0; param.xmax = 1.0;
```

## Extension: Dirichlet bc

We need to assemble on boundary faces. We will need `FEFaceValues` to do this.

```c++
   FEFaceValues<dim> fe_bface_values(fe,
                                     Quadrature<dim-1> (Point<dim>(0.0)),
                                     update_values);
```

In cell loop, distinguish between boundary and other faces

```c++
if(cell->face(0)->at_boundary() && periodic == false)
{
   // assemble boundary flux
   fe_bface_values.reinit(cell, 0);
}
else if(cell->face(1)->at_boundary() && periodic == false)
{
   // assemble boundary flux
   fe_bface_values.reinit(cell, 1);
}
else
{
   // assemble on right face, as we already do
}
```
