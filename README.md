# Codes on finite element method

These codes are prepared as material for teaching finite element method.

* deal.II: Examples using deal.II library in C++
* fenics: Examples using Fenics in Python (fenics is dead)
* dolfinx: rewriting fenics codes in dolfinx
* firedrake: rewriting fenics codes in firedrake, code is on colab
* dg1d: 
  * python: A simple 1-D DG code in Python scalar/euler equations
  * c: A simple 1-D DG code euler equations (very old, not well written)
  * Other directories contain deal.II codes.
* dg2d: 2-D codes written in deal.II

## How to get the code ?

(1) The easiest way is by using git which will download the code in directory `fem`

```shell
git clone https://github.com/cpraveen/fem.git
```

(2) Download a zip file by clicking on the green button called "Code" and then "Download ZIP".

(3) Or use wget

```shell
wget https://github.com/cpraveen/fem/archive/master.zip
```

or curl

```shell
curl -O https://github.com/cpraveen/fem/archive/master.zip
```

and unzip it

```shell
unzip master.zip
```

which creates the directory ```fem-master```.

## Other codes

* [fem50](https://github.com/cpraveen/fem50): matlab codes for 2d BVP
* [juliafem](https://github.com/cpraveen/juliafem): julia codes for 2d BVP
* [deal.II](https://bitbucket.org/cpraveen/deal_ii): more deal.II examples
* [fenics](https://github.com/cpraveen/fenics): more fenics examples

## Additional software

Some of the examples require additional software.

* [Gmsh](http://gmsh.info): used for generating unstructured grids
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/executables): used to visualize solutions 
* [ParaView](https://www.paraview.org): another software to visualize solutions
