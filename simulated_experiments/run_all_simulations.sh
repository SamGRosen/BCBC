sbatch --export=ALL,METHOD='ADAPTIVE_BCBC' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCBC' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCBC_NO_SCALE' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCBC_NO_SCALE_OR_ADAPT' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCEL' run_simulated_trial.sh
sbatch --export=ALL,METHOD='COBRA' run_simulated_trial.sh
