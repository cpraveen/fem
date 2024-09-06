# 2d linear advection equation

> Copied from bitbucket version bd9ec3a7b13e914fd1641f4271dafa3ab691eb7b
> 
> https://bitbucket.org/cpraveen/deal_ii/src/master/dg/2d_scalar_unsteady_unstr/

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

You can save in high order format by setting

```shell
flags.write_higher_order_cells = true
```

but this only works in paraview.

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
