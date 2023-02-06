#!/bin/bash
#SBATCH --time=06:10:00
#SBATCH --account=def-mideon
#SBATCH --cpus-per-task=1
#SBATCH --mem=187G
#SBATCH --mail-type=ALL
#SBATCH --array=1-1100

bash parallel_burn_trunc.sh SIM burn_convex 5 53 $SLURM_ARRAY_TASK_ID burn 7368 0.0000166667 1

cd /scratch/peter114/PopDataProcessing/SIM_VCFs/SIM_burn/
