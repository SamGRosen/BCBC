#!/bin/bash
#SBATCH --output=./slurm-output/%A.out
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

module load R/4.4.0
module load GLPK/5.0
module load OpenBLAS/3.23

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/x86_64-linux-gnu/

pwd

echo $SLURM_JOB_ID
echo $METHOD
echo $DATASET

Rscript /hpc/home/sgr26/BCBC/real_data_experiments/data_experiment.R $SLURM_JOB_ID \
  $METHOD $DATASET
