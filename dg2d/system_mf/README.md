# 2d system of conservation laws on unstructured grids

Matrix-free DG solver for 2d systems of conservation laws using Lagrange basis.
Works on Cartesian and unstructured quadrilateral grids. Gauss-Legendre and
Gauss-Lobatto-Legendre points are used for the basis and quadrature.

Available models:

* `models/euler`: Compressible Euler equations
* `models/linadv`: Linear advection

## Run Euler (isentropic vortex)

```shell
cd fem/dg2d/system_mf
ln -sf ./models/euler/pde.h
ln -sf ./models/euler/isentropic_vortex/problem.h
make release
make
cd models/euler/isentropic_vortex
mpirun -np 4 ../../../main input.prm
```

## Run Linear advection (rotating gaussian)

```shell
cd fem/dg2d/system_mf
ln -sf ./models/linadv/pde.h
ln -sf ./models/linadv/rotate.h problem.h
make release
make
cd models/linadv
mpirun -np 4 ../../main input.prm
```

## Viewing results

```shell
visit -o solution.xdmf
```
