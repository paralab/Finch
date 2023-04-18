#!/bin/bash
#SBATCH --time=0:50:00     # walltime, abbreviated by -t
#SBATCH -o ./outerr/out-job-%j-node-%N.tsv      # name of the stdout, using the job number (%j) and 
                           # the first node (%N)
#SBATCH -e ./outerr/err-job-%j-node-%N.log      # name of the stderr, using job and first node values

#SBATCH --nodes=10         #
#SBATCH --ntasks=512      # maximum total number of mpi tasks across all nodes.
#SBATCH -p normal

#SBATCH --job-name=matvec_bench_channel-weak512_p2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=masado@cs.utah.edu


##
## Note: Run this from the build directory.
##

# load appropriate modules

module load petsc/3.13
module load intel/19.0.5

RUNPROGRAM=/home1/07803/masado/Dendro-KT/build/matvecBenchChannel

PTS_PER_PROC_WEAK=131072

# If there are not enough levels then we may get errors.
MAX_DEPTH=20

ELE_ORDER=2
ibrun -n   1  -o   0 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np001.tsv" &
ibrun -n   2  -o   1 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np002.tsv" &
ibrun -n   4  -o   3 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np004.tsv" &
ibrun -n   8  -o   7 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np008.tsv" &
ibrun -n  16  -o  15 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np016.tsv" &
ibrun -n  32  -o  31 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np032.tsv" &
ibrun -n  64  -o  63 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np064.tsv" &
ibrun -n 128  -o 127 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np128.tsv" &
ibrun -n 256  -o 255 task_affinity $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np256.tsv" &
wait

ibrun -n 512 $RUNPROGRAM $PTS_PER_PROC_WEAK $MAX_DEPTH $ELE_ORDER 4 "simpleDA_p2" > "out-${SLURM_JOB_ID}-np512.tsv"

cat "out-${SLURM_JOB_ID}-"np{001,002,004,008,016,032,064,128,256,512}.tsv | sort -k2 -t'	' -n > "out-${SLURM_JOB_ID}-simpleDA_p2-np001-512.tsv"
rm "out-${SLURM_JOB_ID}-"np{001,002,004,008,016,032,064,128,256,512}.tsv

