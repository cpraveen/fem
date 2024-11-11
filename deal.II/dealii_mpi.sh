# Where do you want to install deal.II, you can change this.
DEAL_II_DIR=$HOME/deal.II

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
   -DEAL_II_ALLOW_AUTODETECTION=OFF \
   -DEAL_II_FORCE_AUTODETECTION=OFF \
   -DEAL_II_ALLOW_BUNDLED=ON \
   -DEAL_II_WITH_VTK=OFF \
   ..
