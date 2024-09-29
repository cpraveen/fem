# 2d system of conservation law

Solves 2d system of conservation laws on Cartesian grids using Legendre basis. The domain need not be a rectangle, but can be anything, but the cells need to be rectangles.

* `models/euler`: Compressible Euler equations

## Run an example

```shell
ln -s ../models/euler/pde.h
ln -s ../models/euler/isentropic_vortex.h problem.h
cmake .
make release
make
./main > log.txt 2>&1 &
tail -f log.txt
visit -o sol*.vtu
```

## Exercise: Maxwells' equations

See Hesthaven, Section 6.5.

## Exercise: Shallow water equations

Implement shallow water equation model and solve radial dam-break problem, see LeVeque, Sec. 21.7.

## Exercise: Linear advection equation

We can also solve scalar conservation law with this code. Implement linear advection equation and solve some IVP.
