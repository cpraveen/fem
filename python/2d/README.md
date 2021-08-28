# 2D Poisson equation using linear FEM on triangles

You need [Gmsh](http://www.gmsh.info) and [meshio](https://github.com/nschloe/meshio) to run this problem.

## Generate the mesh

```bash
gmsh -2 mesh_tri.geo
```

See the mesh

```bash
python mesh_tri.py
```

which looks like this

<p align="center">
<img src="output/mesh.svg">
</p>

The code also shows how to make use of mesh object and how to draw contour plots on a triangular grid.

## Solve the Poisson equation

```bash
python poisson2d.py
```

The computed solution is shown below

<p align="center">
<img src="output/sol.svg">
</p>