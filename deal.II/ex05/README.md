# Solve Poisson equation in 2d with mixed BC

Similar to ex02, same exact solution is used, but with 
  * Dirichlet BC on x=0 and x=1
  * Neumann BC on y=0 and y=1

So we solve the following problem
```
-Laplace(u) = f in (0,1) x (0,1)
u = uexact on x=0 and x=1
u = g      on y=0 and y=1
```
The values for the boundary condition are obtained from the exact solution, and g = d(uexact)/dn
