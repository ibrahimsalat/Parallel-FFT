#!/bin/bash
# =================
# serial_scaling.sh
# =================

#SBATCH --job-name=serial_scaling
#SBATCH --partition=teach_cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:5:00
#SBATCH --mem-per-cpu=15000M

# Load modules required for runtime e.g.
module load languages/intel/2020-u4

cd $SLURM_SUBMIT_DIR
for i in {15..26}
do 
    ./fft $i
done