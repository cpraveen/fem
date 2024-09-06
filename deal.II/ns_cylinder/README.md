# Incompressible flow past a cylinder in channel

Create mesh

```shell
gmsh -2 turek.geo
```

Set some parameters in `turek.prm` file.

Run steady solver

```shell
./main -p turek.prm -steady
```

Run unsteady solver

```shell
./main -p turek.prm -unsteady
```

Run unsteady solver starting from steady solution

```shell
./main -p turek.prm -unsteady -restart
```
