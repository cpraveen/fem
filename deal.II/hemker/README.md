# Convection-diffusion: Hemker problem

This has Galerkin and SUPG.

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
python surface.py
python plot.py
```

Run SUPG

```shell
./demo -supg
python surface.py
python plot.py
```
