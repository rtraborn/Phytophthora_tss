#!/bin/bash

#SBATCH -n 1
#SBATCH -t 0-8:00
#SBATCH -A rraborn
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.out

myDir=/scratch/rraborn/Phytop_tss/scripts/tsr

echo "Launching job"

cd $myDir

./xrunSwf > err

echo "Job complete"
