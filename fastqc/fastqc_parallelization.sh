#!/bin/sh
#SBATCH --time=4:00:00
#SBATCH --array=1-316
#SBATCH --ntasks=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END

sed "${SLURM_ARRAY_TASK_ID}q;d" fastqc.cmdlist | bash