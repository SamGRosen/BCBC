library(tidyverse)
source("simulated_experiments/simulation_helpers.R")

all_checkers <- generate_checkers(2024)

pattern <- "BCBC_NO_SCALE_15517188"

print(paste("Getting paths for pattern", pattern))

paths_to_combine <- list.files(paste0("/cwork/sgr26/", pattern),
                               pattern=pattern, full.names = TRUE)

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
