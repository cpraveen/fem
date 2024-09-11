#!/bin/bash

# Set deal.II version
V=9.6.0

# Where do you want to install deal.II, you can change this.
DEAL_II_DIR=$HOME/deal.II

# For parallel compiling, set this to number of cores in your cpu
# If you have a multicore machine, set a larger value here.
# The compilation will be faster.
NPROC=2

# Needed if you are behind proxy
# If not, then comment/remove following two lines.
#http_proxy=http://192.168.1.1:3128
#https_proxy=$http_proxy

#-------------------------------------------------------------------------------
# DONT CHANGE ANYTHING BELOW THIS LINE UNLESS YOU ARE AN EXPERT
#-------------------------------------------------------------------------------

# Exit on error
set -e

echo "To install with MPI, you need mpicc, mpic++, mpif90"
read -p "Install with MPI ? (y/n/ctr-c) " MPI
if [[ "$MPI" == "y" ]]; then
   if ! command -v  mpicc &> /dev/null
   then
      echo "mpicc not found"
      MPI=n
   fi
   if ! command -v  mpic++ &> /dev/null
   then
      echo "mpic++ not found"
      MPI=n
   fi
   if ! command -v  mpif90 &> /dev/null
   then
      echo "mpif90 not found"
      MPI=n
   fi
   if [[ "$MPI" == "n" ]]; then
      echo "Missing MPI"
      exit
   fi
fi

mkdir -p $DEAL_II_DIR
mkdir -p $DEAL_II_DIR/dealii-build
cd $DEAL_II_DIR/dealii-build

#echo "==> Downloading deal.II cmake script"
#wget -c https://raw.githubusercontent.com/cpraveen/fembook/master/deal.II/dealii.sh

echo "==> Downloading deal.II sources"
wget -c https://github.com/dealii/dealii/releases/download/v${V}/dealii-${V}.tar.gz
echo "==> Extracting deal.II sources"
tar zxvf dealii-${V}.tar.gz > install.log
cd dealii-${V} && rm -rf build && mkdir -p build && cd build
echo "==> Run cmake"

if [[ "$MPI" == "n" ]]; then
# Without MPI
cmake \
   -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
   -DDEAL_II_CXX_FLAGS="-march=native -mtune=native -std=c++17" \
   -DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
   -DDEAL_II_WITH_MPI=OFF \
   -DDEAL_II_WITH_VTK=OFF \
   -DDEAL_II_COMPONENT_EXAMPLES=ON  \
   -DDEAL_II_COMPILE_EXAMPLES=OFF \
   ..
fi

if [[ "$MPI" == "y" ]]; then
# With MPI
cmake \
   -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
   -DDEAL_II_CXX_FLAGS="-march=native -mtune=native -std=c++17" \
   -DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
   -DDEAL_II_WITH_MPI=ON \
   -DCMAKE_C_COMPILER=mpicc  \
   -DCMAKE_CXX_COMPILER=mpic++  \
   -DCMAKE_Fortran_COMPILER=mpif90  \
   -DDEAL_II_COMPONENT_EXAMPLES=ON  \
   -DDEAL_II_COMPILE_EXAMPLES=OFF \
   ..
fi

echo "==> Compiling"
make -j $NPROC
echo "==> Install"
make install

# Download and install documentation
cd $DEAL_II_DIR
echo "==> Download documentation"
wget -c https://github.com/dealii/dealii/releases/download/v${V}/dealii-${V}-offline_documentation.tar.gz
echo "==> Extract documentation"
tar zxvf dealii-${V}-offline_documentation.tar.gz > install.log
#rm dealii-${V}-offline_documentation.tar.gz

echo "==> Installed in $DEAL_II_DIR"
echo "==> Add the following line to your $HOME/.bashrc file if you use bash"
echo "==> or to the startup file for your shell"
echo "    export DEAL_II_DIR=$DEAL_II_DIR"
echo "==> You may delete directory: $DEAL_II_DIR/dealii-build"
