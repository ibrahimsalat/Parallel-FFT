#!/bin/bash
# =================
# mt_scaling.sh
# =================

#SBATCH --job-name=mt_scaling
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --cpus-per-task=1
#SBATCH --time=0:50:00
#SBATCH --mem=0
#SBATCH --exclusive

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

cd $SLURM_SUBMIT_DIR
for i in {15..26}
do 
    ./fft_mt $i 25
done