#!/bin/bash
#SBATCH --time=0:48:00     # walltime, abbreviated by -t
#SBATCH -o ../results/ipdps21/out-job-%j-node-%N.tsv      # name of the stdout, using the job number (%j) and 
                           # the first node (%N)
#SBATCH -e ../results/ipdps21/err-job-%j-node-%N.log      # name of the stderr, using job and first node values

#SBATCH --nodes=293         #
#SBATCH --ntasks=16384      # maximum total number of mpi tasks across all nodes.
#SBATCH -p normal

#SBATCH --job-name=matvec_bench_channel-weak16K_p2
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
ibrun -n 4096  -o    0 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER "weak_p2" > "out-${SLURM_JOB_ID}-weak_p2-np4096.tsv" &
ibrun -n 8192  -o 4096 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER "weak_p2" > "out-${SLURM_JOB_ID}-weak_p2-np8192.tsv" &
wait

ibrun -n 16384 $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER "weak_p2" > "out-${SLURM_JOB_ID}-weak_p2-np16384.tsv"

cat "out-${SLURM_JOB_ID}-weak_p2-"np{4096,8192,16384}.tsv | sort -k2 -t'	' -n > "out-${SLURM_JOB_ID}-weak_p2-np4K-16K.tsv"
rm "out-${SLURM_JOB_ID}-weak_p2-"np{4096,8192,16384}.tsv

