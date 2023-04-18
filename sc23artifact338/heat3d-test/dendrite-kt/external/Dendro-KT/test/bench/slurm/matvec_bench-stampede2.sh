#!/bin/bash
#SBATCH --time=0:20:00     # walltime, abbreviated by -t
#SBATCH -o ../results/out-job-%j-node-%N.tsv      # name of the stdout, using the job number (%j) and 
                           # the first node (%N)
#SBATCH -e ../results/err-job-%j-node-%N.log      # name of the stderr, using job and first node values

#SBATCH --nodes=10         # request 4 nodes, that's the most we'll use
#SBATCH --ntasks=480      # maximum total number of mpi tasks across all nodes.

#SBATCH --job-name=matvec_bench
#SBATCH --mail-type=ALL
#SBATCH --mail-user=masado@cs.utah.edu

# additional information for allocated clusters
#SBATCH --account=TG-DPP130002    # account - abbreviated by -A
#SBATCH --partition=skx-normal  # partition, abbreviated by -p

##
## Note: Run this in the 'results' directory.
##

# load appropriate modules
#module load intel impi

RUNPROGRAM=/home1/03727/tg830270/Research/Dendro-KT/build/matvecBench

# The size of a 4D level==3 regular grid.
PTS_PER_PROC=1024

# The starting level is 3, we need up to two more levels during weak scaling.
# If there are not enough levels then we may get errors.
MAX_DEPTH=$((3+2))


## =============
#Strong scaling.
## =============
ELE_ORDER=1
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $(($PTS_PER_PROC / $NP)) $MAX_DEPTH $ELE_ORDER "strong";
  echo ;
done

ELE_ORDER=2
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $(($PTS_PER_PROC / $NP)) $MAX_DEPTH $ELE_ORDER "strong";
  echo ;
done

ELE_ORDER=4
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $(($PTS_PER_PROC / $NP)) $MAX_DEPTH $ELE_ORDER "strong;"
  echo ;
done


## =============
#Weak scaling.
## =============
ELE_ORDER=1
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $PTS_PER_PROC $MAX_DEPTH $ELE_ORDER "weak";
  echo ;
done

ELE_ORDER=2
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $PTS_PER_PROC $MAX_DEPTH $ELE_ORDER "weak";
  echo ;
done

ELE_ORDER=4
for NP in 1 16 256 ; do
  ibrun -np $NP $RUNPROGRAM $PTS_PER_PROC $MAX_DEPTH $ELE_ORDER "weak";
  echo ;
done
