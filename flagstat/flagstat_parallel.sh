#!/bin/sh
#SBATCH --time=12:00:00
#SBATCH --array=1-379
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END
#SBATCH --mail-user=dwight.newton@mail.utoronto.ca



sed "${SLURM_ARRAY_TASK_ID}q;d" flagstat.cmdlist | bash


