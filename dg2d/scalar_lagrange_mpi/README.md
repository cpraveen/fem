# 2d linear advection equation

> Copied from bitbucket version bd9ec3a7b13e914fd1641f4271dafa3ab691eb7b
> 
> https://bitbucket.org/cpraveen/deal_ii/src/master/dg/2d_scalar_unsteady_unstr/

This can work on both Cartesian and quad meshes. It can run in parallel using parallel vectors from Petsc or Trilinos. Deal.II must have been compiled with Petsc/Trilinos for this to work.

In this code, we use `MeshWorker::loop` for assembly unlike in some other codes where we used `MeshWorker::mesh_loop`.

## Run the code

Set `grid` in `run` function.

For `grid = 2`, generate grid

```shell
gmsh -2 annulus.geo
```

Compile

```shell
make release
make
```

Run

```shell
mpirun -np 4 ./dg
```

See solution

```shell
visit -o solution.visit
```

or

```shell
paraview solution.pvd
```

You can save the solution in high order format by setting

```shell
flags.write_higher_order_cells = true
```

but this only works in paraview.

With Visit, you can see the mesh partitions by plotting `Add -> Subset -> domains`.

## Exercise

Write a function

```c++
template <int dim>
void refine_grid(Triangulation<dim> &triangulation, double r0, double r1)
{
   // Fill this function
}
```

which refines all cells whose vertices lie between radial distance `r0 <= r <= r1`. After the grid is set up, call this function two times

```c++
refine_grid(triangulation, 0.25, 0.75);
refine_grid(triangulation, 0.40, 0.60);
```

See `step-1` in the deal.II tutorial. A solution is given in the file `refine_grid.cc` but do not look at it until you have attempted yourself.

## TODO

The code solves linear advection equation and this is hard-coded. Make this more general by extracting the PDE part to an include file as we have done in other examples.
