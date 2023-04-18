# TalyFEM's Dependencies

## PETSc (3.6.x or 3.7.x)

```bash
# Download PETSc
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.6.4.tar.gz
tar xvf petsc-lite-3.6.4.tar.gz
rm petsc-lite-3.6.4.tar.gz

# Configure PETSc (will also download and compile BLAS/LAPACK, ParMETIS, and some other requirements). Requires a Fortran compiler for one of the PETSc modules (fblaslapack).
cd petsc-3.6.4
./configure --with-clanguage=cxx --download-scalapack --download-mpich --download-parmetis --download-fblaslapack --download-blacs --download-metis --download-mumps

# Set PETSC_DIR to get automatically set on login via .bashrc
# (this is how CMake knows where to find PETSc)
echo "export PETSC_DIR=`pwd`" >> $HOME/.bashrc
echo "export PETSC_ARCH=arch-linux2-cxx-debug" >> $HOME/.bashrc

# Reload .bashrc (so PETSC_DIR is set for the next step)
source ~/.bashrc

# compile
make
```

## libconfig

```bash
# Download and extract
wget http://www.hyperrealm.com/packages/libconfig-1.4.10.tar.gz
tar xvf libconfig-1.4.10.tar.gz
rm libconfig-1.4.10.tar.gz

# Compile and write files to libconfig-1.4.10/install
cd libconfig-1.4.10
./configure --prefix=`pwd`/install
make
make install

# Permanently set $LIBCONFIG_DIR environment variable
echo "export LIBCONFIG_DIR=`pwd`/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

## HDF5 (Optional)

```bash
# Download it
wget https://support.hdfgroup.org/ftp/HDF5/current18/src/hdf5-1.8.19.tar

# Extract it and remove the archive (we're done with it)
tar xvf hdf5-1.8.19.tar
rm hdf5-1.8.19.tar

# Compile HDF5
# This will install the binary to the `pwd`/install directory.
cd hdf5-1.8.19
CC=mpicc ./configure --enable-parallel --prefix=`pwd`/install
make
make install
```

# TalyFEM

## Desktop

```bash
sudo apt-get install git cmake3 libconfig-dev build-essential gcc g++ gfortran
# install PETSc separately (see start of guide)

# Clone the TalyFEM repository
git clone https://bitbucket.org/baskargroup/taly_fem
cd taly_fem

# Build
mkdir build
cd build
cmake ..
make -j8

# Test (but skip the 1,500 SSHT tests)
ctest -j4 -E ssht
```

## CyEnce

```bash
# There is a shared PETSc pre-compiled on CyEnce (in optimized mode)
echo "export PETSC_DIR=/work/baskargroup/shared/petsc-3.6.4-intel17" >> ~/.bashrc
echo "export PETSC_ARCH=arch-linux2-cxx-opt" >> ~/.bashrc

# Same for libconfig
echo "export LIBCONFIG_DIR=/work/baskargroup/shared/libconfig-1.5/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc

echo "module load intel/17.0.1" >> ~/.bashrc

# the version of CMake on CyEnce by default is ancient (2.x)
echo "module load cmake/3.5.0-rc2" >> ~/.bashrc

# Apply our changes to the current environment
source ~/.bashrc

# finally, try and build TalyFEM
git clone https://bitbucket.org/baskargroup/taly_fem
cd taly_fem
mkdir build
cd build
cmake ..
make
```

## Condo

```bash
# There is a shared PETSc pre-compiled on Condo
echo "export PETSC_DIR=/work/baskargroup/shared/petsc-3.6.4-working" >> ~/.bashrc
echo "export PETSC_ARCH=arch-linux2-cxx-debug" >> ~/.bashrc

# Same for libconfig
echo "export LIBCONFIG_DIR=/work/baskargroup/shared/libconfig-1.5-intel17/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc

echo "module load cmake/3.5.0-rc2" >> ~/.bashrc
echo "module load intel/17.0.1" >> ~/.bashrc

source ~/.bashrc

# finally, try and build TalyFEM
git clone https://bitbucket.org/baskargroup/taly_fem
cd taly_fem
mkdir build
cd build
cmake ..
make
```

## Comet

Unfortunately, you will need to compile PETSc yourself on Comet. The available module was compiled using the Intel 13 compiler, which is too old to compile TalyFEM with (missing C++11 support), and the newer compiler's runtime is not backwards-compatible with programs compiled with the Intel 13 compiler. If you don't, you'll get strange linker errors about undefined references to `intel_sse2_*` symbols.

```bash
# The default intel module is 2013, which is too old
echo "module swap intel intel/2016.3.210" >> ~/.bashrc
echo "module load cmake" >> ~/.bashrc

# Apply changes
source ~/.bashrc

# install PETSc, make sure you do ./configure --with-cc=mpicc --with-cxx=mpicxx
# install libconfig, make sure you use ./configure CC=mpicc CXX=mpicxx --prefix=`pwd`/install

# finally, try and build
mkdir build
cd build
cmake ..
make
```
