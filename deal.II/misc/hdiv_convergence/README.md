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
    0    64   144 1.403e-02    - 3.552e-01    - 
    1   256   544 3.516e-03 2.00 1.780e-01 1.00 
    2  1024  2112 8.797e-04 2.00 8.902e-02 1.00 
    3  4096  8320 2.200e-04 2.00 4.452e-02 1.00 
    4 16384 33024 5.499e-05 2.00 2.226e-02 1.00 

==============================================================
 Convergence Results for RT(1) Space in 2D  Cartesian grid
==============================================================

cycle cells  dofs     L2_error      Hdiv_error   
    0    64    544 3.476e-04    - 1.802e-02    - 
    1   256   2112 4.351e-05 3.00 4.511e-03 2.00 
    2  1024   8320 5.441e-06 3.00 1.128e-03 2.00 
    3  4096  33024 6.802e-07 3.00 2.821e-04 2.00 
    4 16384 131584 8.502e-08 3.00 7.053e-05 2.00 

==============================================================
 Convergence Results for RT(0) Space in 2D  Quadrilateral grid
==============================================================

cycle cells dofs     L2_error      Hdiv_error   
    0    64   144 1.842e-02    - 3.758e-01    - 
    1   256   544 6.929e-03 1.41 2.229e-01 0.75 
    2  1024  2112 3.279e-03 1.08 1.667e-01 0.42 
    3  4096  8320 1.674e-03 0.97 1.519e-01 0.13 
    4 16384 33024 8.419e-04 0.99 1.481e-01 0.04 

==============================================================
 Convergence Results for RT(1) Space in 2D  Quadrilateral grid
==============================================================

cycle cells  dofs     L2_error      Hdiv_error   
    0    64    544 5.779e-04    - 2.049e-02    - 
    1   256   2112 1.399e-04 2.05 7.487e-03 1.45 
    2  1024   8320 3.200e-05 2.13 2.931e-03 1.35 
    3  4096  33024 7.797e-06 2.04 1.365e-03 1.10 
    4 16384 131584 1.962e-06 1.99 6.726e-04 1.02 
```
