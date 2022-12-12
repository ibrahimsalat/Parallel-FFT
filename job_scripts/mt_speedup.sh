#!/bin/bash
# =================
# mt_speedup.sh
# =================

#SBATCH --job-name=mt_speedup
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
for i in {1..28}
do 
    ./fft_mt 26 $i
done