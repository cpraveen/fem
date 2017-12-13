# Set deal.II version
V=8.5.1

# Where do you want to install deal.II
DEAL_II_DIR=$HOME/deal.II

# Needed if you are behind proxy
# If not, then comment/remove following two lines.
http_proxy=http://192.168.1.1:3128
https_proxy=$http_proxy

mkdir -p $DEAL_II_DIR

cd $HOME
mkdir -p dealii-build
cd dealii-build
echo "==> Download deal.II cmake script"
wget -c https://raw.githubusercontent.com/cpraveen/fembook/master/deal.II/dealii.sh

echo "==> Downloading deal.II sources"
wget -c https://github.com/dealii/dealii/releases/download/v${V}/dealii-${V}.tar.gz
echo "==> Extracting deal.II sources"
tar zxvf dealii-${V}.tar.gz > install.log
cd dealii-${V}
mkdir -p build
cd build
echo "==> Run cmake"
sh $HOME/dealii-build/dealii.sh
echo "==> Compiling"
# If you have a multicore machine, set a larger value here. 
# The compilation will be faster.
make -j 2
echo "==> Install"
make install

cd $DEAL_II_DIR
echo "==> Download documentation"
wget -c https://github.com/dealii/dealii/releases/download/v${V}/dealii-${V}-offline_documentation.tar.gz
echo "==> Extract documentation"
tar zxvf dealii-${V}-offline_documentation.tar.gz > install.log
rm dealii-${V}-offline_documentation.tar.gz

echo "==> Installed in $DEAL_II_DIR"
echo "==> Add this to your $HOME/.bashrc file"
echo "==> export DEAL_II_DIR=$DEAL_II_DIR"
echo "==> You may delete directory: $HOME/dealii-build"
