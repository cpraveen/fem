# Error convergence under uniform refinement

This code solves Poisson equation on unit square on succesively refined grids.

```shell
cmake .
make
./demo
```

This produces a latex file `error.tex`

```shell
pdflatex error.tex
```

Open and see `error.pdf` for a table of error versus mesh size and corresponding convergence rate.

The solution error as a function of x,y is saved into different files as the grid is refined. You can open the error-##.vtk files in visit and see the error distribution in the domain.

## Discontinuous dirichlet bc

Solve the problem

```text
   -Laplace(u) = 0  in (-0.5, +0.5) x (0, 1)
```

with dirichlet bc which is obtained from the exact solution

```text
   u = (1/pi) * atan(y/x)
```

On the bottom side `y=0`, the dirichlet bc is discontinuous at the origin. Change the code like this. The exact solution:

```c++
template <>
double ExactSolution<2>::value (const Point<2> &p,
                                const unsigned int /*component*/) const
{
   const double x = p[0];
   const double y = p[1];
   return (1.0/M_PI) * atan2(y, x);
}

template<>
Tensor<1,2> ExactSolution<2>::gradient (const Point<2>   &p,
                                        const unsigned int) const
{
   Tensor<1,2> values;
   const double x = p[0];
   const double y = p[1];
   const double r2 = std::pow(x,2) + std::pow(y,2);
   values[0] = -y/(M_PI * r2);
   values[1] =  x/(M_PI * r2);
   return values;
}
```

The right hand side:

```c++
template <>
double RightHandSide<2>::value (const Point<2> &p,
                                const unsigned int /*component*/) const
{
   return 0.0;
}
```

The mesh

```c++
   Point<dim> p0(-0.5, 0.0);
   Point<dim> p1(+0.5, 1.0);
   GridGenerator::hyper_rectangle (triangulation, p0, p1)
```

Now we see that the solution converges at first order only, while derivative does not converge at all. The error is mostly around the origin where the solution is less regular.
