# 1-d linear scalar PDE using DG

This uses nodal Lagrange basis either with Gauss-Legendre or Gauss-Lobatto-Legendre nodes; the same nodes are also used for all quadrature.

To solve constant linear advection equation

```shell
cp ../1d_scalar_legendre/linadv/*.h .
```

To solve Burgers equation

```shell
cp ../1d_scalar_legendre/burgers/*.h .
```

Then, see the readme file in `1d_scalar_legendre` directory.

Select GL or GLL basis in `input.prm` file

```shell
set basis = gl | gll
```
