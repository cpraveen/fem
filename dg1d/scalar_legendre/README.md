# 1-d scalar PDE using DG

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
./main input.prm
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
FunctionParser<1> initial_condition("sin(2*pi*x)","pi=3.141592653589793");
FunctionParser<1> exact_solution("sin(2*pi*(x-t))","pi=3.141592653589793");
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
Functions::SymbolicFunction<1> exact_solution("sin(2*pi*(x-t))");
exact_solution.set_time(param.final_time);
param.xmin = 0.0; param.xmax = 1.0;
```

## Exercise: Dirichlet bc

Solve the following problem

```text
u_t + u_x = 0 for x in (0,1)
IC   : u(x,0) = sin(2*pi*x)
BC   : u(0,t) = sin(-2*pi*t)
Exact: u(x,t) = sin(2*pi*(x-t))
```

Do not apply periodicity to the mesh.

We need to assemble flux on boundary faces and we will use upwind flux.

In cell loop, distinguish between boundary and other faces

```c++
if(cell->face(0)->at_boundary() && param->periodic == false)
{
   // assemble left boundary flux
}

if(cell->face(1)->at_boundary() && param->periodic == false)
{
   // assemble right boundary flux
}
else
{
   // assemble on right face, as we already do
}
```

## Exercise: Classical RK4 scheme

Implement the classical four stage, fourth order RK scheme.
