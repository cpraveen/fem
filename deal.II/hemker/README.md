# Convection-diffusion: Hemker problem

See

> Hemker, A singularly perturbed model problem for numerical computation, J. Comp. App. Math., vol. 76, pp. 277-285, 1996. https://doi.org/10.1016/S0377-0427(96)00113-6

The code has Galerkin and SUPG methods.

Generate mesh

```shell
gmsh -2 hemker.geo
```

Compile

```shell
make release
make
```

Run Galerkin

```shell
./demo -eps 1.0e-3
python plot.py
python surface.py (rotate with mouse)
```

There are lot of oscillations in the solution.

Run Galerkin with refined mesh

```shell
./demo -eps 1.0e-3 -nrefine 1
python plot.py
python surface.py (rotate with mouse)
```

The oscillations are still present.

Run SUPG

```shell
./demo -eps 1.0e-3 -supg 1
python plot.py
python surface.py (rotate with mouse)
```

There are no noticeable oscillations in the solution.

## Refining the grid in gmsh

To refine within Gmsh,

```shell
gmsh hemker.geo
Mesh -> 2D
Mesh -> Refine by splitting  # Repeat this for more refinement
Mesh -> Save
```
