# Codes on finite element method
These codes are prepared as material for teaching finite element method. Most of the codes are based on deal.II library

## deal.II
The examples are based on deal.II finite element library. You can find detailed instructions on the deal.II website. To run most of the examples given here, it is enough to compile deal.II in serial mode. You will need a C++ compiler and cmake.

Suppose we install deal.II in `$HOME/deal.II` directory. Add following line in your `$HOME/.bashrc` file
```
export DEAL_II_DIR=$HOME/deal.II
```
To compile deal.II, follow these steps. A sample `dealii.sh` script is included in this repository.
```
cd $HOME
wget https://github.com/dealii/dealii/releases/download/v8.4.2/dealii-8.4.2.tar.gz
tar zxvf dealii-8.4.2.tar.gz
cd dealii-8.4.2
mkdir build
cd build
sh /path/to/dealii.sh
make -j4
make install
```
Now you can delete the directory where you compiled and also the source file
```
cd $HOME
rm -rf dealii-8.4.2
rm dealii-8.4.2.tar.gz
```
Also, download and install the offline documentation by following these steps.
```
cd $HOME/deal.II
wget https://github.com/dealii/dealii/releases/download/v8.4.2/dealii-8.4.2-offline_documentation.tar.gz
tar zxvf dealii-8.4.2-offline_documentation.tar.gz
rm dealii-8.4.2-offline_documentation.tar.gz
```
Now you can open `$HOME/deal.II/doc/index.html` in your web browser and view the documentation.
