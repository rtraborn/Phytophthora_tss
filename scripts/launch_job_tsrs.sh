#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-8:00
#SBATCH -p cmecpu1
#SBATCH -q cmeqos
#SBATCH -A rraborn
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.out

myDir=/scratch/rraborn/Phytop_tss/scripts/tsr

module load r/4.0.2

echo "Launching job"

cd $myDir

./xrunSwf > err

echo "Job complete"
