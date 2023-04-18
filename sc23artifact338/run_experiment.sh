#!/bin/bash

# This script assumes Julia has been installed and can
# be run with the command "julia"

# At least 55 CPU cores and 8 GPUs must be allocated.

# On Frontera, the following modules should be installed:
module load intel/19.1.1 impi/19.0.9 petsc/3.16

# Libconfig is needed if recompiling
# export LIBCONFIG_DIR=$(pwd)/libconfig-1.7.2/install
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBCONFIG_DIR/lib

# Set up the packages
julia setup.jl;

# Run the Finch script to generate the code
julia heat3d-dendrite.jl;
# The code is now in the heat3d folder

# To use, copy the contents of heat3d/ into the 
# folder with dendrite-kt, make a build folder,
# use cmake with this command:
#    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc
# then make with:
#    make
# and be sure to copy the config.txt file into build/

# To simplify the process for this test, use a 
# code that is already in place and includes
# an extra output specifically for this test.
cd heat3d-test

# If running on Frontera with the modules loaded, 
# this is already compiled.

# remove previous results
rm solvetime.txt

# Run (estimated 1.5 hours)
for i in 512 256 128 64 32 16 8 4 2 1;
do
    echo "Running with " $i " processes";
    ibrun -n $i ./heat3d
done

# The times are recorded in solvetime.txt
cd ..

# Process the results
julia process.jl;

# The generated figure is in results.png
# The VTK files are in heat3d-test/result_00001/