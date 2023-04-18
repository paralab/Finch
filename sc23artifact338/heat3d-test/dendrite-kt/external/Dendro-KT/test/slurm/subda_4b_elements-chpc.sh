#!/bin/bash
#SBATCH --time=6:00:00     # walltime, abbreviated by -t
#SBATCH -o out-job-%j-node-%N.out      # name of the stdout, using the job number (%j) and 
                                       # the first node (%N)
#SBATCH -e err-job-%j-node-%N.err      # name of the stderr, using job and first node values

#SBATCH --nodes=6          # request 6 node
#SBATCH --ntasks=96        # maximum total number of mpi tasks across all nodes.

#SBATCH --job-name=DendroKT-4b
#SBATCH --mail-type=ALL
#SBATCH --mail-user=masado@cs.utah.edu

# additional information for allocated clusters
#SBATCH --account=soc-kp    # account - abbreviated by -A
#SBATCH --partition=soc-kp  # partition, abbreviated by -p


# load appropriate modules
module load intel impi

RUNPROGRAM=./../build/tstBillionsOfElements

mpirun -np $SLURM_NTASKS ./$RUNPROGRAM
