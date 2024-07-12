#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --array=1-2

module load R/4.4.0
module load GLPK/5.0
module load OpenBLAS/3.23

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/

pwd

Rscript /hpc/home/sgr26/BCBC/test.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_JOB_ID $SLURM_JOB_ID $SLURM_ARRAY_TASK_ID
