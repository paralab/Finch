This document will show how to obtain and compile TalyFEM and its dependencies.

# Before You Start

You will need:

0. A Linux system ([Ubuntu](https://en.wikipedia.org/wiki/Ubuntu_(operating_system)) is recommended) and a basic understanding of how to use the terminal.

1. [git](https://en.wikipedia.org/wiki/Git), a popular version control system. On Ubuntu, you can install git by running `sudo apt-get install git` from the terminal.

    Most HPC clusters already have git installed by default. The version of git you use is typically not important.

2. [cmake](https://cmake.org/), version 3.5 or above. CMake is a popular cross-platform build system for C++ projects. Like git, it can be installed on Ubuntu by running `sudo apt-get install cmake` in a terminal. Use `cmake --version` to see what version of CMake is currently installed.

    Many HPC clusters have a version of CMake available by default, but it is often too out of date to use with TalyFEM. On HPC systems, you can usually use `module avail cmake` and `module load cmake/[version here]` to try loading a newer version. If that fails, you will have to install it yourself - download a binary distribution from the CMake website.

3. C, C++, and Fortran compilers. Your C++ compiler must support the C++11 standard. On Ubuntu, the easiest way to get all of these is `sudo apt-get install build-essential gcc g++ gfortran`.

    Most HPC clusters will have these available by default, but the C++ compiler is often out of date. You can look at `module avail` to see what is available. (Don't try to install a compiler yourself on an HPC system - you should be able to find a module.)

# Required Dependencies

TalyFEM has two required dependencies:

1. [PETSc](https://www.mcs.anl.gov/petsc/), a mature C library of routines for solving partial differential equations in parallel at large scales. We mostly use it for solving the sparse `Ax=B` systems that the finite element method generates in parallel. PETSc uses [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) to operate in parallel.

2. [libconfig](http://www.hyperrealm.com/oss_libconfig.shtml), a mature C/C++ library for parsing config files. TalyFEM has built-in utilities for loading simulation parameters from a file, and this library is what we use to parse that file. At the moment, these config files are tightly coupled with TalyFEM, so this dependency is required.


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
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.6.4.tar.gz
tar xvf petsc-lite-3.6.4.tar.gz
rm petsc-lite-3.6.4.tar.gz  # delete the archive, we're done with it

# Configure PETSc, downloading and compiling dependencies along the way.
cd petsc-3.6.4
./configure --download-mpich --download-metis --download-parmetis --download-scalapack --download-fblaslapack --download-mumps

# Write commands to set the PETSC_DIR and PETSC_ARCH environmnt variables to our .bashrc,
# which is automatically executed on login. These environment variables are what TalyFEM
# uses to find the PETSc installation.
echo "export PETSC_DIR=`pwd`" >> ~/.bashrc
echo "export PETSC_ARCH=arch-linux2-c-debug" >> ~/.bashrc

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
wget http://www.hyperrealm.com/packages/libconfig-1.4.10.tar.gz
tar xvf libconfig-1.4.10.tar.gz
rm libconfig-1.4.10.tar.gz

# Compile and copy output files to libconfig-1.4.10/install
cd libconfig-1.4.10
./configure --prefix=`pwd`/install
make -j4  # compile with 4 threads
make install

# Permanently set $LIBCONFIG_DIR environment variable, which is what TalyFEM uses
# to find your libconfig install. Also set LD_LIBRARY_PATH to include the
# libconfig lib directory to prevent dynamic linking errors with libconfig++.so.
echo "export LIBCONFIG_DIR=`pwd`/install" >> ~/.bashrc
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$LIBCONFIG_DIR/lib" >> ~/.bashrc
source ~/.bashrc
```

**NOTE:** On some HPC clusters (notably Condo when using the Intel compiler), the `make` step gives a linker error. This is the libconfig example program failing to link with the Intel runtime. This is okay - the libconfig library itself compiles just fine. Just run `make install` and double check that `install/lib` contains some `*.a` files.

# Compiling TalyFEM

TalyFEM uses [CMake](https://cmake.org/) as its build system. CMake will first look in the standard system directories for dependencies (i.e. `/usr`), and then check the directories specified with the `$PETSC_DIR`, `$PETSC_ARCH`, and `$LIBCONFIG_DIR` environment variables if they are set. You can specify the compiler with `-DCMAKE_CXX_COMPILER=mpicxx`.

**ACTUALLY IMPORTANT NOTE:** CMake caches some values between runs. If you are having CMake configuration issues, clean the build directory before you re-run CMake (`rm -rf build/*`).

```bash
# Clone the TalyFEM repository
git clone https://bitbucket.org/baskargroup/taly_fem
cd taly_fem

# Build everything inside the "build" subdirectory
mkdir build
cd build
cmake ..  # the ".." tells CMake where to find our CMakeLists.txt (up a directory)
make -j8  # compile with 8 threads in parallel
```

See [docs/CMAKE.md](CMAKE.md) for more info on how to use CMake.

## Enabling HDF5 I/O

TalyFEM optionally supports a custom [HDF5](https://en.wikipedia.org/wiki/Hierarchical_Data_Format)-based binary input/output format, which can be visualized using ParaView (via XDMF files). To enable HDF5 support, the `ENABLE_HDF5` CMake variable must be set to `YES`, and CMake must be able to find a working installation of the HDF5 C library. Example:

```bash
export HDF5_ROOT=/path/to/hdf5/install
cmake .. -DENABLE_HDF5=YES
```

(`export HDF5_ROOT=...` is only necessary if your HDF5 installation is in a non-standard location)

(the CMake FindHDF5 module checks the HDF5_ROOT *environment variable*, not the CMake variable!)

## Using PTSCOTCH for Domain Decomposition

TalyFEM can optionally use the [SCOTCH](https://www.labri.fr/perso/pelegrin/scotch/) library for parallel domain decomposition instead of ParMETIS. It is enabled similarly to HDF5:

```bash
cmake .. -DENABLE_SCOTCH=YES -DSCOTCH_DIR=/path/to/scotch
```

(`-DSCOTCH_DIR=...` is only necessary if your SCOTCH installation is in a non-standard location)

**NOTE:** The SCOTCH domain decomposition code requires the SCOTCH ParMETIS shim to be compiled with a prefix (since we use both ParMETIS and PTSCOTCH simultaneously). You must compile PTSCOTCH with the `-DSCOTCH_METIS_PREFIX` option (this is a SCOTCH option, not a TalyFEM option!).

**NOTE:** The size of SCOTCH indices must match the size of `PetscInt`. `PetscInt` is 32-bit by default, but SCOTCH indices are 64-bit by default. So, unless you are using 64-bit PETSc indices, you must compile SCOTCH with `-DIDXSIZE32` (this is a SCOTCH option, not a TalyFEM option!).

**NOTE:** ParMETIS is currently still used to convert the FEM mesh to a dual graph for partitioning, so you will still need ParMETIS installed.

## Enabling MESQUITE

TalyFEM includes an optional [MESQUITE](https://trilinos.org/oldsite/packages/mesquite/) mesh adapter (`TalyMesqMesh`). This allows one to use MESQUITE's mesh improvements algorithms directly on a TalyFEM GRID.
Note that there is no automatic use of MESQUITE. You will need to manually call the MESQUITE mesh optimization routines yourself. See `tests/mesquite/src/main.cpp` for an example.

It is enabled by:

```bash
cmake .. -DENABLE_MESQUITE=YES -DMESQUITE_DIR=/path/to/mesquite
```

(`-DMESQUITE_DIR` is only encessary if your MESQUITE installation is in a non-standard location. You can alternatively set the `$MESQUITE_DIR` environment variable (e.g. in your bashrc) to avoid specifying it when running CMake every time. In either case, `$MESQUITE_DIR/include/Mesquite.hpp` and `$MESQUITE_DIR/lib/libmesquite.la` should exist.)

Note that the mesh adapter currently only works for non-domain-decomposed meshes, and in the parallel-but-not-DD case it is possible that non-deterministic MESQUITE optimizers could give different meshes on each process. If you intend to use the MESQUITE adapter at scale, you should investigate these issues first.

## Using 64-bit PETSc Indices (problems with >= 2^31 unknowns)

TalyFEM will work out of the box with a PETSc that has been compiled with the `--with-64-bit-indices` option. No extra configuration is necessary. 64-bit indices are only necessary for very large problems where the size of the matrix/vector system is larger than a 32-bit signed integer can represent (`2^31 - 1`).

**NOTE:** The MUMPS solver is not supported with 64-bit indices, so make sure you remove `--download-mumps` from your PETSc configure options.

## Generating the Doxygen Documentation

If [Doxygen](http://www.stack.nl/~dimitri/doxygen/) was installed when you configured TalyFEM, you can run `make docs` to generate HTML documentation for TalyFEM and the tutorial projects. Open `build/docs/html/index.html` in a web browser to take a look. The docs for the tutorials will be in `build/tutorials/[tutorial_name]/docs/html/index.html`.

## Testing TalyFEM

TalyFEM also includes a number of system and unit tests (`tutorials/*` and `tests/*`). These are hooked up to CTest, the CMake testing framework. You can run the tests in serial with `make test`, or in parallel using `ctest` with the `-j` option:

```bash
ctest -j4
```

This will spawn 4 threads to chew through all available tests. There are a number of other useful CTest options, including options to re-run only failed tests, resume an interrupted test run, and include or exclude tests by name. Read the [CTest manual](https://cmake.org/cmake/help/latest/manual/ctest.1.html) for more information.

**NOTE:** Some tests will test TalyFEM's parallel features, but this information has not been conveyed to CTest. These kinds of tests may spawn up to 4 processes, so it's a good idea to be conservative with your `-j` parameter.

**NOTE:** The steady-state heat tutorial exhaustively tests boundary condition combinations, and that takes a long time. You can run with `ctest -j4 -E ssht` to exclude these (1,500 or so) tests.

## Static Analysis with clang-tidy

If you are using CMake 3.7 or above, TalyFEM can automatically run [clang-tidy](http://clang.llvm.org/extra/clang-tidy/) on source code during compilation. Clang-tidy uses LLVM (the backend for the clang compiler) to parse your source code and provide additional code correctness and style warnings. The current set of checks is quite conservative, but may be expanded in the future.

If the `clang-tidy` executable is found in the normal system paths during configuration, TalyFEM will automatically use it. On Ubuntu, you can install clang-tidy with `sudo apt-get install clang-tidy` (or `clang-tidy-3.9` on older versions of Ubuntu). Note that using clang-tidy can make library compilation take nearly 1.5 times as long (but the warnings it provides could save you some debugging).

# Creating a TalyFEM Project

First, you'll need a CMakeLists.txt. You can copy one from one of the existing tutorials/tests, or use this one:

```cmake
# This is a required line that specifies the minimum version of CMake
# our CMakeLists.txt requires.
cmake_minimum_required(VERSION 3.5.0 FATAL_ERROR)

# This tells CMake what the name of our project is (optional).
project(HelloWorld)

# This declares an executable named 'hello' and lists the .cpp files that need
# to be compiled to make the executable. Different files are separated by
# whitespace (not commas).
add_executable(hello
  # list source files to compile
  src/main.cpp
)

# Next, we add the ./include directory to the include search path for the program.
# This is what lets us "#include" any header files we put in the "include/" folder.
# ${CMAKE_CURRENT_SOURCE_DIR} is a CMake variable that corresponds to the
# directory containing CMakeLists.txt.
target_include_directories(hello PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Now we tell CMake to find TalyFEM.
# The "REQUIRED" means if it can't be found stop and print an error.
find_package(talyfem REQUIRED)

# Finally, we link our excecutable defined above with TalyFEM. This will
# also link our program with TalyFEM's dependencies (PETSc, libconfig).
target_link_libraries(hello talyfem)
```

Update the source files list accordingly as you add more .cpp files. (At least on Windows, you may also need to add .h files for them to show up inside your IDE's project tree.) To build, run:

```bash
mkdir build
cd build
cmake .. -Dtalyfem_DIR=/path/to/taly/build  # tell CMake what directory contains talyfemConfig.cmake
make
```

You *do not* need to re-run CMake every time you modify your project. The `Makefile` that CMake generates will automatically re-run CMake when necessary.

**IMPORTANT:** CMake caches a lot of things inside its build directory. If you make drastic changes to your CMake command-line arguments after the initial run, you may need to clean your build directory with `rm -rf build/*` before running CMake again. You likely need to do this if you want to change your C++ compiler (`-DCMAKE_CXX_COMPILER=...`).

**REALLY IMPORTANT:** The generated CMake build files are not relative! You cannot move a build directory to a new project directory and expect it to compile the code inside the new project - the generated build files are hard-coded to use absolute paths to the old directory. You must regenerate your CMake directory with the `cmake` command (as shown above).

**LESS IMPORTANT:** You should probably not commit anything in your build directory to version control. To prevent accidents, you can create a [.gitignore](https://www.atlassian.com/git/tutorials/gitignore) file to tell git to ignore the build directory. A .gitignore file is just a text file with each line containing patterns for files that should be ignored by git (e.g. `build/`).
