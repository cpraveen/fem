# 1-D DG code for scalar conservation laws

Solves linear advection and Burger's equation with periodic boundary condition. The basis functions are Legendre polynomials and quadrature is performed using Gauss-Legendre points. To get help, type
```
python dg.py -h
```

## Linear advection equation

Try a smooth initial condition
```
python dg.py -pde linear -ic sine -degree 1 -ncell 50
```
After any multiple of time period, the solution is equal to initial condition and we can compute the L2 error norm
```
python dg.py -pde linear -ic sine -degree 1 -ncell 50 -compute_error yes
```
Now try with a discontinuous initial condition
```
python dg.py -pde linear -ic hat -degree 1 -ncell 52
```
To suppress the oscillations, we will use a limiter
```
python dg.py -pde linear -ic hat -degree 1 -ncell 52 -limit yes
```

## Inviscid Burgers equation

Try a smooth initial condition
```
python dg.py -pde burger -ic sine -degree 1 -ncell 50
```
Discontinuous initial condition
```
python dg.py -pde burger -ic hat -degree 1 -ncell 52
```
Discontinuous initial condition with limiter
```
python dg.py -pde burger -ic hat -degree 1 -ncell 52 -limit yes
```
