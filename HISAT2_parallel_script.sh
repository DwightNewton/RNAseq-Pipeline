#!/bin/sh
#SBATCH --time=96:00:00
#SBATCH --array=1-80
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END

sed "${SLURM_ARRAY_TASK_ID}q;d" hisat_human.cmdlist | bash