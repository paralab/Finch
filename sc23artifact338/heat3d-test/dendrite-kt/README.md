This document will show how to obtain and compile dendrite-kt and its dependencies.

**Cloning the repository**
```bash
git clone https://bitbucket.org/baskargroup/dendrite-kt.git
git submodule update --init
```

# Required Dependencies

dendrite-kt has these required dependencies:

1. A C++ compiler with C++11 support is required. We aim for compatibility with gcc and icc (Intel C++ compiler).

2. [PETSc >= 3.7, upto 3.13.4 tested](https://www.mcs.anl.gov/petsc/), a mature C library of routines for solving partial 
differential equations in parallel at large scales. We mostly use it for solving the sparse `Ax=B` 
systems that the finite element method generates in parallel. 
PETSc uses [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) to operate in parallel.

3. [libconfig](http://hyperrealm.github.io/libconfig/), a mature C/C++ library for parsing 
config files. Dendrite-kt has built-in utilities for loading simulation parameters from a file, 
and this library is what we use to parse that file. 
At the moment, these config files are tightly coupled with Dendrite-kt, so this dependency is required.

4. [BLAS](http://www.netlib.org/blas/)

5. [LAPACK](http://www.netlib.org/lapack/)

6. [LAPACKE](https://www.netlib.org/lapack/lapacke.html)

## Compiling PETSc

PETSc has its own build system which is a combination of Python and CMake. To build PETSc, you'll run `./configure` with some options, and then run `make` to build the PETSc library.

PETSc is built on top of a lot of other libraries. Thankfully, the PETSc build system has options to automatically download and install these libraries if desired. This guide recommends the following options for compiling PETSc:

```bash
./configure --download-mpich --download-metis --download-parmetis --download-scalapack --download-fblaslapack --download-mumps
```

Explanations for each argument:

* `--download-mpich`: download [MPICH](https://www.mpich.org/). MPI is a popular standard in scientific computing for writing parallel programs, and MPICH is one implementation of the standard. MPICH is maintained by the same people who maintain PETSc.

* `--download-metis` and `--download-parmetis`: download [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview), a graph partitioning library, as well as its parallel version, ParMETIS. This is used by TalyFEM to do domain decomposition.

* `--download-scalapack` and `--download-fblaslapack`: these are some basic linear algebra routines that PETSc uses to implement some of its built-in solvers. Compiling fblaslapack requires a working Fortran compiler.

* `--download-mumps`: download the MUMPS direct linear solver. This is an alternative to the built-in direct solver, and it works in parallel. It's optional.

How to install PETSc in the terminal:

```bash
# Download and extract PETSc
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.13.4.tar.gz
tar xvf petsc-lite-3.13.4.tar.gz
rm petsc-lite-3.13.4.tar.gz  # delete the archive, we're done with it

# Configure PETSc, downloading and compiling dependencies along the way.
cd petsc-3.13.4
./configure --download-mpich --download-metis --download-parmetis --download-scalapack --download-fblaslapack --download-mumps

# Write commands to set the PETSC_DIR and PETSC_ARCH environmnt variables to our .bashrc,
# which is automatically executed on login. These environment variables are what TalyFEM
# uses to find the PETSc installation.
echo "export PETSC_DIR=`pwd`" >> ~/.bashrc
echo "export PETSC_ARCH=arch-linux-c-debug" >> ~/.bashrc
source ~/.bashrc

# To add mpicc and mpicxx to the default PATH
echo "export PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH" >> ~/.bashrc
# Reload .bashrc now so PETSC_DIR and PETSC_ARCH are set.
source ~/.bashrc

# Compile PETSc
# (this will be automatically multi-threaded, so we don't need -j)
make
```

**NOTE:** Some HPC systems already have an MPI implementation available (e.g. Intel MPI or Cray MPI). On these systems you *may* get better performance by using this system default MPI instead of MPICH. You can do this by removing the `--download-mpich` option and specifying to use the system-provided MPI compiler wrapper as the PETSc compiler (e.g. with `--with-cc=mpicc` and `--with-cxx=mpicxx`). PETSc should automatically pick up the MPI directory from the MPI compiler. If it doesn't, you can try to specify it manually with `--with-mpi-dir=/some/path/here`.

**NOTE:** PETSc passes along compiler flags to automatically add its libraries to the [rpath](https://en.wikipedia.org/wiki/Rpath) of your compiled project executable. This is why you don't need to modify `LD_LIBRARY_PATH` in this step.

**NOTE:** For production runs, you can specify `--with-debugging=0` to turn off PETSc's internal debugging features. This may give a performance boost (at the cost of less-comprehendable error messages). Do not do any program development using a version of PETSc compiled this way.

## Compiling libconfig

Libconfig's build system is [Autotools](https://www.gnu.org/software/automake/manual/html_node/Autotools-Introduction.html), which means you'll need to run `./configure` and then `make` to build.

This guide recommends passing ```--prefix=`pwd`/install``` to `./configure`, which will cause `make install` to copy the output files to `[your_libconfig_dir]/install` instead of `/usr`. This way your libconfig install lives completely inside your libconfig folder. This is necessary if you are working on a system where you don't have admin privileges (i.e. an HPC cluster).

```bash
# Download and extract
wget http://hyperrealm.github.io/libconfig/dist/libconfig-1.7.2.tar.gz
tar xvf libconfig-1.7.2.tar.gz
rm libconfig-1.7.2.tar.gz

# Compile and copy output files to libconfig-1.7.2/install
cd libconfig-1.7.2
./configure --prefix=`pwd`/install
make -j8  # compile with 8 threads
make install

# Permanently set $LIBCONFIG_DIR environment variable, which is what TalyFEM uses
# to find your libconfig install. Also set LD_LIBRARY_PATH to include the
# libconfig lib directory to prevent dynamic linking errors with libconfig++.so.
echo "export LIBCONFIG_DIR=`pwd`/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

**NOTE:** On some HPC clusters (notably Condo when using the Intel compiler), the `make` step gives a linker error. This is the libconfig example program failing to link with the Intel runtime. This is okay - the libconfig library itself compiles just fine. Just run `make install` and double check that `install/lib` contains some `*.a` files.



Compilation
============
```bash
mkdir build;
cd build;
```
### For 2D  Computation
```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_2D=Yes -DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
make
```

### For 3D  Computation
```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_3D=Yes -DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
make
```

### For 4D  Computation
```bash
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DENABLE_4D=Yes -DCMAKE_C_COMPILEER=mpicc -DCMAKE_CXX_COMPILER=mpicxx
make
```

### For IBM Computation (WIP)
Additional CMake Flags that is required:
```bash
-DIBM = Yes
```

### For Optimized Matrix Assembly (WIP) 
Additional CMake Flags that is required:
```bash
-DTENSOR=Yes
```

## Compile using CLion
We recommand using CLion as IDE. To compile dendrite-kt under CLion environment:

1. add a toolchain in ```File | Settings | Build, Execution, Deployment | Toolchains```, 
set ```C Compiler=$(which mpicc)``` and ```C++ Compiler=$(which mpicxx)```. 

2. In ```File | Settings | Build, Execution, Deployment | CMake```, change ```Toolchain``` to the one we just created, 
add environment variable ```PETSC_DIR``` and ```PETSC_ARCH``` in ```environment```, if CMake is complaining about not 
finding ```Libconfig++```, add ```LIBCONFIG_DIR``` and ```LD_LIBRARY_PATH``` also.

3. Add additional ```CMAKE_FLAGS``` in ```CMAKE Options```

Development in CLion will offer multiple choice of ```Toolchain``` and multiple ```Profile``` that makes it easier 
for changing build options or switch between different executables.

## Compile under ubuntu-20 focal

1. Change gcc version from 9 to 7

2. Install python2 if not available

Testing Dendrite-kt
============
Dendrite-kt also includes a number of regression tests, they ensure that the library is correctly installed and any
changes to the library does not break existing features.
To run the regression test, [python3](https://www.python.org/about/) needs to be installed.
### Running individual regression test
make sure the 2D and 3D executables are compiled 
```bash
cd examples/SSHT
python3 ssht_regression_test.py --runner=$(which mpirun)
```
or 
```bash
cd examples
sh regression_test.sh SSHT
```
The second approach is preferred as it will automatically compile the executables.

### Running all tests
This will run all the regression tests, took roughly 5~10 minutes.
```bash
cd examples
sh regression_test.sh
```

Contributors
============
* Masado Ishii
* Kumar Saurabh
* Hari Sundar
* Baskar Ganapathysubramanian
