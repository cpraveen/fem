# Codes on finite element method
These codes are prepared as material for teaching finite element method.

* deal.II: Examples using deal.II library in C++
* fenics: Examples using Fenics in Python
* dg1d: A simple 1-D DG code in Python for scalar conservation law

## How to get the code ?

(1) The easiest way is by using git which will download the code in directory ```fembook```
```
git clone https://github.com/cpraveen/fembook.git
```
(2) Download a zip file by clicking on the green button called "Clone or download" and then "Download ZIP".

(3) Or use wget
```
wget https://github.com/cpraveen/fembook/archive/master.zip
```
and unzip it
```
unzip master.zip
```
which creates the directory ```fembook-master```.

## Additional software
Some of the examples require additional software.

* [Gmsh](http://gmsh.info): used for generating unstructured grids
* [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/executables): used to visualize solutions 
* [ParaView](https://www.paraview.org): another software to visualize solutions
