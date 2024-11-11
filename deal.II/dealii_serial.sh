# You can change this location to something else if you want.
DEAL_II_DIR=$HOME/deal.II

# Without MPI
cmake \
   -DCMAKE_INSTALL_PREFIX=$DEAL_II_DIR \
   -DDEAL_II_CXX_FLAGS="-march=native -mtune=native -std=c++17" \
   -DDEAL_II_CXX_FLAGS_RELEASE="-O3" \
   -DDEAL_II_WITH_MPI=OFF \
   -DDEAL_II_COMPONENT_EXAMPLES=ON  \
   -DDEAL_II_COMPILE_EXAMPLES=OFF \
   -DEAL_II_ALLOW_AUTODETECTION=OFF \
   -DEAL_II_FORCE_AUTODETECTION=OFF \
   -DEAL_II_ALLOW_BUNDLED=ON \
   -DEAL_II_WITH_VTK=OFF \
   ..
