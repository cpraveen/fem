# 2d system of conservation law on Cartesian grids

Solves 2d system of conservation laws on Cartesian grids using Legendre basis. The domain need not be a rectangle, but can be anything, but the cells need to be rectangles with sides aligned with x, y axes.

* `models/euler`: Compressible Euler equations

## Run an example

In a terminal, go to `fem/dg2d/system_legendre`

```shell
ln -s ../models/euler/pde.h
ln -s ../models/euler/isentropic_vortex/problem.h
cmake .
make release
make
```

In another terminal, go to `fem/dg2d/models/euler/isentropic_vortex`

```shell
../../../system_legendre/main input.prm > log.txt 2>&1 &
tail -f log.txt
visit -o sol*.vtu
```

## Exercise: Linear advection equation

We can also solve scalar conservation law with this code. Implement linear advection equation and solve some IVP, see `scalar_legendre` code. This is done in `models/linadv` but try to do it yourself before seeing that solution.

## Exercise: Maxwells' equations

See Hesthaven, Section 6.5.

## Exercise: Shallow water equations

Implement shallow water equation model and solve radial dam-break problem, see LeVeque, Sec. 21.7.

## Exercise: Implement multi-threading

We use MeshWorker to assemble the DG residual which automatically uses multiple threads depending on the number of cores. You can use `top` command and examine the CPU usage, and you will see more than 100% usage if multiple threads are being used. You can control the number of threads by setting the shell variable, e.g., to use 6 threads

```shell
export DEAL_II_NUM_THREADS=6
```

Time the program with single and multiple threads

```shell
time ./main > log.txt 2>&1 &
```

You can use MeshWorker for applying TVD limiter and for computing cell average. The TVD limiter is expensive so it might benefit from multi-threading. You can try to use the default scratch data and copy data objects

https://dealii.org/current/doxygen/deal.II/classMeshWorker_1_1ScratchData.html  
https://dealii.org/current/doxygen/deal.II/structMeshWorker_1_1CopyData.html

You have to implement a cell worker function that takes one `cell` and applies the limiter on that cell.
