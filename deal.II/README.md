# FEM examples using deal.II

The codes in directory `deal.II` are based on [deal.II](http://www.dealii.org) library. This is a C++ library for FEM in 1/2/3 dimensions. The examples here do not have detailed explanations and will not be sufficient to learn deal.II by themselves. Please study the superbly documented examples on the deal.II [website](https://www.dealii.org/developer/doxygen/deal.II/Tutorial.html). This is also a great resource to learn C++. You can also consult the following books for more on C++:

1. [B. Stroustrup: Tour of C++](http://www.stroustrup.com/tour2.html)
1. P. Gottschling: Discovering Modern C++

## Examples

1. ex01: 1-d BVP, Dirichlet bc
1. ex02: 2-d Poisson equation on unit square, Dirichlet bc, 2d/3d
   * Quad mesh
   * Triangle mesh
1. ex03: 2-d Poisson equation on unit square, Dirichlet bc, assembly using MatrixCreator
   * Quad mesh
   * Triangle mesh
1. ex04: 2-d Poisson equation on unit square, Dirichlet bc, convergence under uniform grid refinement
1. ex05: 2-d Poisson on unit square, mixed bc
1. ex06: 2-d Laplace equation in non-convex (Gamma-shaped) domain, convergence under uniform refinement
1. ex07: 2-d Laplace equation in non-convex (Gamma-shaped) domain, convergence under adaptive refinement
1. ex08: 2-d Laplace equation with discontinuous Dirichlet bc, uniform and adaptive refinement
1. ex09: 2-d Poisson with discontinuous coefficient. Same as [step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html) in deal.II tutorials, modified to use gmsh grid
1. ex10: Poisson equation in annulus, curved boundaries
1. ex21: 2-d unsteady heat equation
1. ex31: 2-d linear advection equation using Galerkin method
1. ns_cylinder: Incompressible NS for flow over cylinder

## deal.II Examples

1. [step-6](https://dealii.org/developer/doxygen/deal.II/step_6.html): Poisson equation with discontinuous coefficients, grid adaptation
1. [step-10](https://dealii.org/developer/doxygen/deal.II/step_10.html): Higher order mappings, compute value of pi
1. [step-11](https://dealii.org/developer/doxygen/deal.II/step_11.html): Poisson equation with Neumann boundary conditions
1. [step-12](https://dealii.org/developer/doxygen/deal.II/step_12.html): Discontinuous Galerkin methods for steady linear advection problem
1. [step-20](https://dealii.org/developer/doxygen/deal.II/step_20.html): Poisson/Darcy equation using mixed method; Raviart-Thomas/DG
1. [step-22](https://dealii.org/developer/doxygen/deal.II/step_22.html): Stokes flow, Taylor-Hood
1. [step-55](https://dealii.org/developer/doxygen/deal.II/step_55.html): Stokes flow, Taylor-Hood, parallel
1. [step-23](https://dealii.org/developer/doxygen/deal.II/step_23.html): Second order wave equation solved as first order system

## Advanced topics

1. [step-37](https://dealii.org/developer/doxygen/deal.II/step_37.html): Matrix-free approach to Poisson equation

## Installing deal.II

The examples are based on deal.II finite element library. You can find detailed installation instructions on the [deal.II website](http://www.dealii.org/developer/readme.html). To run most of the examples given here, it is enough to compile deal.II in serial mode. You will need a C++ compiler and cmake, which can be installed using the package manager in your operating system. Download the latest stable version of deal.II [here](https://github.com/dealii/dealii/releases). The instructions below are written assuming `v9.6.0` but change this to the actual version you are using.

### Prepare your system

You will need to have C and C++ compilers installed. In Debian/Ubuntu, you can install like this

```shell
sudo apt update
sudo apt install build-essential
```

We need `cmake` to generate makefiles and we also use `wget` below to download some files; also install fortran compiler and gnuplot for visualization,

```shell
sudo apt install gfortran gnuplot cmake wget
```

If you want to run parallel codes, or use Petsc/Trilinos, you need an MPI library. Check if you already have MPI

```shell
which mpirun
```

If not present, then install openmpi (You can also use mpich)

```shell
sudo apt install openmpi-bin openmpi-common openmpi-doc
```

We will also need some visualization software which you should install from the following websites.

* [VisIt](https://visit-dav.github.io/visit-website/index.html)
* [Paraview](https://www.paraview.org)

For grid generation, we need Gmsh, get it from [here](https://www.gmsh.info).

### Automated installation

A bash script that can automate the installation of deal.II is also included. Do the following steps

```shell
cd $HOME
wget https://raw.githubusercontent.com/cpraveen/fem/master/deal.II/dealii_basic.sh
# or use curl if you dont have wget
curl -O https://raw.githubusercontent.com/cpraveen/fem/master/deal.II/dealii_basic.sh
bash ./dealii_basic.sh
```

This will install deal.II into the directory `$HOME/deal.II`; you can change the location of install directory by editing the bash script.

### Manual download and install

Suppose we want to install deal.II in `$HOME/deal.II` directory. Add following line in your `$HOME/.bashrc` file

```shell
export DEAL_II_DIR=$HOME/deal.II
```

You may have to source the file

```shell
source $HOME/.bashrc
```

or open a new terminal window.

To compile deal.II, follow these steps. A sample `dealii.sh` script is included in this repository.

```shell
cd $HOME
wget https://github.com/dealii/dealii/releases/download/v9.6.0/dealii-9.6.0.tar.gz
tar zxvf dealii-9.6.0.tar.gz
cd dealii-9.6.0
wget https://raw.githubusercontent.com/cpraveen/fembook/master/deal.II/dealii.sh
mkdir -p build && cd build
bash ../dealii.sh
```

This creates Makefile etc., and the contents of the `build`  directory look something like this

```shell
$ ls
bin		CMakeFiles	     detailed.log  lib		 source
bundled		cmake_install.cmake  doc	   Makefile	 summary.log
cmake		contrib		     examples	   revision.log  tests
CMakeCache.txt	CTestTestfile.cmake  include	   share
```

Now compile

```shell
make -j4      # -jN,  N = number of CPU cores on your computer
make install
```

It is good to keep the build files in case you want to re-compile later. Otherwise, you can delete the directory where you compiled and also the source file

```shell
cd $HOME
rm -rf dealii-9.6.0
rm dealii-9.6.0.tar.gz
```

Also, download and install the offline documentation by following these steps.

```shell
cd $DEAL_II_DIR
wget https://github.com/dealii/dealii/releases/download/v9.6.0/dealii-9.6.0-offline_documentation.tar.gz
tar zxvf dealii-9.6.0-offline_documentation.tar.gz
rm dealii-9.6.0-offline_documentation.tar.gz
```

Now you can open `$DEAL_II_DIR/doc/index.html` in your web browser and view the documentation.

## Using docker on Linux

Docker images are available

https://hub.docker.com/r/dealii/dealii

Install

```shell
docker pull dealii/dealii
```

Run as

```shell
cd fem
docker run -ti --name dealii -v $(pwd):/home/dealii/shared \
                             -w /home/dealii/shared  \
                             dealii/dealii:latest
```

This creates a container with name `dealii` and you can also see this in dashboard. The directory `fem` is now visible inside this container. If you exit the terminal, the container will stop running. To join the container again, first start the container and then attach

```shell
docker start dealii
docker attach dealii
```

To attach a second terminal to the same running container, do

```shell
docker exec -it dealii /bin/bash
```

You can delete the container using the dashboard or on terminal

```shell
docker remove dealii
```

To reduce your typing work, create an alias in your `$HOME/.bashrc` file

```shell
alias docker_dealii="docker run -ti --name dealii \
                                -v $(pwd):/home/dealii/shared \
                                -w /home/dealii/shared  \
                                dealii/dealii:latest"
```

Now you can start the container by typing in a new terminal

```shell
docker_dealii
```

## Install deal.II using apt

On Debian, Ubuntu and other `apt`-based systems, you may be able to install `deal.II` using `apt` packages, see

https://github.com/dealii/dealii/wiki/Debian-and-Ubuntu 

However, you may not get the latest version of deal.II; in that case, you can install all the dependencies using `apt` and compile `deal.II` yourself.
