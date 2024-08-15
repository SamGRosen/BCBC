library(tidyverse)
source("simulated_experiments/simulation_helpers.R")


set.seed(2024)
all_checkers <- list()

trials <- 20
noise_levels <- c(3)
extra_dims <- c(0, 50, 100, 250, 500, 1000)

i <- 1
for(trial in 1:trials) {
  for (noise_level in noise_levels) {
    for (extra_dim in extra_dims) {
      all_checkers[[i]] <- gen_checkerboard(
        100,
        100,
        5,
        5,
        noise = noise_level,
        cluster_spread = 10,
        p_extra = extra_dim,
        prob_empty = 0
      )
      all_checkers[[i]]$noise_level = noise_level
      all_checkers[[i]]$extra_dim = extra_dim
      all_checkers[[i]]$trial = trial
      i <- i + 1
    }
  }
}

pattern <- "BCBC_NO_SCALE_OR_ADAPT_12940400"

print(paste("Getting paths for pattern", pattern))

paths_to_combine <- list.files("/cwork/sgr26", pattern=pattern, full.names = TRUE)

run_ids <- as.integer(str_match(paths_to_combine, "_(\\d{1,3})\\.RDS")[, 2])

paths_in_order <- paths_to_combine[order(run_ids)]

print(paste("Reading", length(paths_in_order), "RDS files"))

all_rds <- paths_in_order |> map(readRDS)

save_path <- paste0("/cwork/sgr26/", pattern, ".RDS")

print(paste("Saving to", save_path))

saveRDS(all_rds, save_path)

csv_path <- paste0("~/BCBC/simulated_experiments/rds_files/", pattern, ".csv")

print(paste("Saving csv to", csv_path))

as_csv <- get_all_checker_results(all_checkers, all_rds, bcbc_extractor)

write.csv(as_csv, csv_path, row.names=FALSE)
