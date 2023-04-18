#!/bin/bash
#SBATCH --time=0:06:00     # walltime, abbreviated by -t
#SBATCH -o ../results/ipdps21/out-job-%j-node-%N.tsv      # name of the stdout, using the job number (%j) and 
                           # the first node (%N)
#SBATCH -e ../results/ipdps21/err-job-%j-node-%N.log      # name of the stderr, using job and first node values

#SBATCH --nodes=37         #
#SBATCH --ntasks=2048      # maximum total number of mpi tasks across all nodes.
#SBATCH -p development

#SBATCH --job-name=matvec_bench_channel-weak2048_p2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=masado@cs.utah.edu


##
## Note: Run this from the build directory.
##

# load appropriate modules

module load petsc/3.13
module load intel/19.0.5

RUNPROGRAM=/home1/07803/masado/Dendro-KT/build/matvecBenchChannel

PTS_PER_PROC_WEAK=1024

# If there are not enough levels then we may get errors.
MAX_DEPTH=15

ELE_ORDER=2
ibrun -n 1024 $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER "weak_p2" > "out-${SLURM_JOB_ID}-weak_p2-np1024-2048.tsv"
ibrun -n 2048 $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER "weak_p2" >> "out-${SLURM_JOB_ID}-weak_p2-np1024-2048.tsv"
