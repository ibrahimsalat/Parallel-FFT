#!/bin/bash
# =================
# mpi32script.sh
# =================

#SBATCH --job-name=mpi__job
#SBATCH --partition=teach_cpu
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=0:50:00
#SBATCH --mem=0

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

cd $SLURM_SUBMIT_DIR

echo "64 proceses with 28 stages"
srun --mpi=pmi2 ./mpi.exe 26
