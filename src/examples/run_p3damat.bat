#!/bin/bash
#
#SBATCH -M notchpeak
#SBATCH -J amattest
#SBATCH --partition=soc-np
#SBATCH --account=soc-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -t 0:20:00

module load julia cmake/3.7.2
module load intel/2019.5.281  impi/2019.5.281 petsc/3.12.0

# #######################################################
# Modify for your setup
# #######################################################
# Enter your path to aMat and Finch
AMAT=~/aMat
FINCH=~/Finch

# The example script name and generated directory name
EXAMPLE=example-poisson3d-amat.jl
GENDIR=p3damat

# The number of processes/partitions to use
NP=32
# #######################################################

# First remove any previously generated version
rm -r $GENDIR
rm -r $AMAT/$GENDIR

# Run the julia example which will generate the code and data files
mpiexec -n $NP julia $FINCH/src/examples/$EXAMPLE

# move the generated directory to the amat dir
# Note: it was generated wherever pwd is
mv $GENDIR $AMAT/

# Set things up in amat
cd $AMAT
rm CMakeLists.txt
rm -r build
mkdir build

cp $GENDIR/CMakeLists.txt ./
cp $GENDIR/MeshData* build/
cp $GENDIR/GeoData* build/
cp $GENDIR/RefelData build/
cd build

# Make it
cmake ..
cmake --build .

# Run it
mpirun -np $NP ./$GENDIR 1 0 0 out

