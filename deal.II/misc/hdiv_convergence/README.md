# Interpolation error for RT(N) space

## Cartesian grid

```math
\| u - u_h\|_{L^2} = O(h^{N+1}), \qquad
\| div(u) - div(u_h)\|_{L^2} = O(h^{N+1})
```

We observe superconvergence of divergence errors.

## Quadrilateral grid

```math
\| u - u_h\|_{L^2} = O(h^{N+1}), \qquad
\| div(u) - div(u_h)\|_{L^2} = O(h^{N})
```

For $N=0$, the L2 error is optimal but we do not observe any convergence in divergence. See

> Pavel B. Bochev and Denis Ridzal, Rehabilitation of the lowest-order Raviart–Thomas element on quadrilateral grids, SIAM J. Numer. Anal. Vol. 47, No. 1, pp. 487–507. https://doi.org/10.1137/070704265

```text
==============================================================
 Convergence Results for RT(0) Space in 2D  Cartesian grid
==============================================================

cycle cells dofs     L2_error      Hdiv_error   
    0    64   144 8.018e-02    - 5.012e-01    - 
    1   256   544 4.008e-02 1.00 2.515e-01 0.99 
    2  1024  2112 2.004e-02 1.00 1.259e-01 1.00 
    3  4096  8320 1.002e-02 1.00 6.295e-02 1.00 
    4 16384 33024 5.010e-03 1.00 3.148e-02 1.00 

==============================================================
 Convergence Results for RT(1) Space in 2D  Cartesian grid
==============================================================

cycle cells  dofs     L2_error      Hdiv_error   
    0    64    544 4.065e-03    - 2.548e-02    - 
    1   256   2112 1.016e-03 2.00 6.380e-03 2.00 
    2  1024   8320 2.540e-04 2.00 1.596e-03 2.00 
    3  4096  33024 6.350e-05 2.00 3.990e-04 2.00 
    4 16384 131584 1.587e-05 2.00 9.974e-05 2.00 

==============================================================
 Convergence Results for RT(0) Space in 2D  Quadrilateral grid
==============================================================

cycle cells dofs     L2_error      Hdiv_error   
    0    64   144 8.191e-02    - 5.196e-01    - 
    1   256   544 4.069e-02 1.01 2.837e-01 0.87 
    2  1024  2112 2.034e-02 1.00 1.869e-01 0.60 
    3  4096  8320 1.018e-02 1.00 1.568e-01 0.25 
    4 16384 33024 5.089e-03 1.00 1.506e-01 0.06 

==============================================================
 Convergence Results for RT(1) Space in 2D  Quadrilateral grid
==============================================================

cycle cells  dofs     L2_error      Hdiv_error   
    0    64    544 4.157e-03    - 2.845e-02    - 
    1   256   2112 1.042e-03 2.00 1.065e-02 1.42 
    2  1024   8320 2.597e-04 2.00 3.876e-03 1.46 
    3  4096  33024 6.512e-05 2.00 1.908e-03 1.02 
    4 16384 131584 1.627e-05 2.00 9.498e-04 1.01 
```
