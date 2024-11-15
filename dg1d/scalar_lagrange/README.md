# 1-d scalar PDE using DG

This uses nodal Lagrange basis either with Gauss-Legendre or Gauss-Lobatto-Legendre nodes; the same nodes are also used for all quadrature.

The code is very similar in most parts to the `scalar_legendre` code, you can see the differences like this

```shell
code -n -d ../scalar_legendre/dg.h ./dg.h
```

To solve constant linear advection equation

```shell
cp ../scalar_legendre/linadv/*.h .
```

To solve Burgers equation

```shell
cp ../scalar_legendre/burgers/*.h .
```

Then, see the readme file in `scalar_legendre` directory.

Select GL or GLL basis in `input.prm` file

```shell
set basis = gl | gll
```

## TODO: Implement strong form DG
