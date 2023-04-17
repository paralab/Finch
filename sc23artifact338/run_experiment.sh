#!/bin/bash

# This script assumes MPI and Julia are available
# and at least 55 CPU cores and 8 GPUs have been allocated.

# Set up the packages
julia setup.jl;

# remove previous results
rm solvetime.txt

# Run the CPU-only version (estimated 2 hours)
for i in 55 40 20 10 5 2 1;
do
    echo "Running with " $i " processes";
    mpiexec -n $i julia bte2d-cpu.jl;
done

# Run the GPU-accelerated version (estimated 30 min)
for i in 8 4 2 1;
do
    echo "Running with " $i " processes/GPUs";
    mpiexec -n $i julia bte2d-gpu.jl;
done

# Process the results
julia process.jl;

# The generated figure is in results.png