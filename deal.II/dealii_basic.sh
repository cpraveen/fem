# Set deal.II version
V=9.2.0

# Where do you want to install deal.II
DEAL_II_DIR=$HOME/deal.II

# For parallel compiling, set this to number of cores in your cpu
NPROC=2

# Needed if you are behind proxy
# If not, then comment/remove following two lines.
#http_proxy=http://192.168.1.1:3128
#https_proxy=$http_proxy

# DONT CHANGE ANYTHING BELOW THIS LINE

mkdir -p $DEAL_II_DIR
mkdir -p $DEAL_II_DIR/dealii-build
cd $DEAL_II_DIR/dealii-build

echo "==> Download deal.II cmake script"
wget -c https://raw.githubusercontent.com/cpraveen/fembook/master/deal.II/dealii.sh

echo "==> Downloading deal.II sources"
wget -c https://github.com/dealii/dealii/releases/download/v${V}/dealii-${V}.tar.gz
echo "==> Extracting deal.II sources"
tar zxvf dealii-${V}.tar.gz > install.log
cd dealii-${V}
rm -rf build 
mkdir -p build
cd build
echo "==> Run cmake"
sh $DEAL_II_DIR/dealii-build/dealii.sh
echo "==> Compiling"
# If you have a multicore machine, set a larger value here. 
# The compilation will be faster.
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
echo "==> Add the following line to your $HOME/.bashrc file"
echo "    export DEAL_II_DIR=$DEAL_II_DIR"
echo "==> You may delete directory: $DEAL_II_DIR/dealii-build"
