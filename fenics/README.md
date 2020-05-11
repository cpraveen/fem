# Examples using Fenics

[Fenics](https://fenicsproject.org) can be installed by following the instructions [here](https://fenicsproject.org/download). For more learning resources, see
* [Fenics tutorial](https://fenicsproject.org/tutorial/)
* [Fenics book](https://fenicsproject.org/book/)

Jupyter notebook version of these examples are available [here](https://drive.google.com/drive/folders/1uonODhFK6gjouBDXNYXHnBAz-nDzqlJL?usp=sharing). I try to maintain identical code in both, but if you find some discrepancy between the two, please inform me.

**THESE CODES ARE TESTED WITH FENICS VERSION 2019.1.0, THEY WILL NOT WORK WITH DOLFINX**

* step-01
  * Laplace equation with constant data
  * Dirichlet BC

* step-02
  * Laplace equation with variable data
  * Dirichlet BC
  * Expression 
  * saving to file and Paraview

* step-03
  * Laplace equation with variable data
  * Dirichlet BC
  * Grid refinement and error convergence

* step-04
  * Laplace equation
  * Dirichlet and Neuman BC
  * Defining different boundary parts using class
  * Mark boundary facets
  * Define new measure
  * Setting up bc list

* step-05
  * Laplace equation in Gamma-shaped domain
  * Use Gmsh to generate mesh
  * Convert to xml
  * Study error convergence by grid refinement

* step-06
  * Laplace equation with discontinuous dirichlet BC

* step-07
  * Laplace equation on unit square
  * pure Neumann bc

* step-08
  * Laplace equation in unit square
  * Dirichlet BC
  * assemble matrices and solve
  * project gradient
  * interpolate
  * max error using numpy

* step-09
  * Laplace equation (same as in step-06)
  * grid adaptation using a-posteriori error indicator

* step-10
  * Laplace equation with discontinuous coefficient
  * grid adaptation using a-posteriori error indicator

* step-11
  * Heat equation in 1-d
  * demo1.py:
    * With time dependent source term
    * Homogeneous Dirichlet bc
  * demo2.py:
    * With time dependent Dirichlet bc

* step-12
  * Steady scalar advection in unit square: a . grad(u) = 0
  * Galerkin, SUPG and DG scheme
