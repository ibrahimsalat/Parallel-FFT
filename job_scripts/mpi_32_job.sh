#!/bin/bash
# =================
# mpi32script.sh
# =================

#SBATCH --job-name=mpi_32_job
#SBATCH --partition=teach_cpu
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=0:50:00
#SBATCH --mem=0

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

cd $SLURM_SUBMIT_DIR
echo "32 proceses with 26 stages"
srun --mpi=pmi2 ./mpi.exe 26
