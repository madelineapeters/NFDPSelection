#!/bin/bash
#SBATCH --time=05:10:00
#SBATCH --account=def-mideon
#SBATCH --cpus-per-task=1
#SBATCH --mem=187G
#SBATCH --mail-type=ALL
#SBATCH --array=2-1100

bash parallel_burn.sh SIM dummy_constant 5 50 $SLURM_ARRAY_TASK_ID burn  0.0000166667 1

cd /scratch/peter114/PopDataProcessing/SIM_VCFs/SIM_burn/
