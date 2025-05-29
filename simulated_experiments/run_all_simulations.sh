sbatch --export=ALL,METHOD='BCBC_NO_SCALE' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCBC_NO_SCALE_OR_ADAPT' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCBC_ADAPT_APPROX' run_simulated_trial.sh
sbatch --export=ALL,METHOD='BCEL' run_simulated_trial.sh
sbatch --export=ALL,METHOD='COBRA' run_simulated_trial.sh

sbatch --export=ALL,METHOD='BCBC_NO_SCALE' run_simulated_trial_equal_space.sh
sbatch --export=ALL,METHOD='BCBC_NO_SCALE_OR_ADAPT' run_simulated_trial_equal_space.sh
sbatch --export=ALL,METHOD='BCBC_ADAPT_APPROX' run_simulated_trial_equal_space.sh
sbatch --export=ALL,METHOD='BCEL' run_simulated_trial_equal_space.sh
sbatch --export=ALL,METHOD='COBRA' run_simulated_trial_equal_space.sh
