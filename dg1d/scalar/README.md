# 1-D DG code for scalar conservation laws

This is a Python code which solves linear advection and Burger's equation with periodic boundary condition and requires Python3, Numpy and Matplotlib.. The basis functions are Legendre polynomials and quadrature is performed using Gauss-Legendre points. To get help, type

```
python dg.py -h
```

## Linear advection equation

Try a smooth initial condition

```
python dg.py -pde linear -ic sin2pi -degree 1 -ncell 50
```

After any multiple of time period, the solution is equal to initial condition and we can compute the L2 error norm

```
python dg.py -pde linear -ic sin2pi -degree 1 -ncell 50 -compute_error yes
```

Now try with a discontinuous initial condition

```
python dg.py -pde linear -ic hat -degree 1 -ncell 52
```

To suppress the oscillations, we will use a limiter

```
python dg.py -pde linear -ic hat -degree 1 -ncell 52 -limit yes
```

The TVD limiter applied to linear advection will cause clipping of smooth extrema

```
python dg.py -pde linear -ic sin2pi -degree 1 -ncell 50 -limit yes
```

## Inviscid Burgers equation

Try a smooth initial condition

```
python dg.py -pde burger -ic sin2pi -degree 1 -ncell 50
```

Discontinuous initial condition

```
python dg.py -pde burger -ic hat -degree 1 -ncell 52
```

Discontinuous initial condition with limiter

```
python dg.py -pde burger -ic hat -degree 1 -ncell 52 -limit yes
```

## Exercises

* Use interpolation to set initial condition. You can use inverse of Vandermonde matrix to convert from nodal values to modal values.
* Try the test cases with a central flux and comment on the quality of the solutions.
* Solve the linear advection equation on a sequence of meshes and compute the rate at which the error decreases with respect to the cell size dx.
* Solve variable coefficient linear advection equation ```u_t + (a(x) u)_x = 0``` with ```a(x) = 1 + (1-x^2)^5``` on [-1,+1] and compare with solution in Hesthaven and Warburton, fig. 5.1
* Rewrite the Python code in any other language of your choice like Fortran, C or C++.
