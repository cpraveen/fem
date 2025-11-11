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
./demo
python surface.py (rotate with mouse)
python plot.py
```

Run SUPG

```shell
./demo -supg
python surface.py (rotate with mouse)
python plot.py
```
