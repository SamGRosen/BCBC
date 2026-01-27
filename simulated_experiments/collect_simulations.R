library(tidyverse)

source("simulated_experiments/simulation_helpers.R")

args <- commandArgs(trailingOnly = TRUE)

path = args[1]

if(grepl("SNR", pattern) > 0) {
  all_checkers <- generate_checkers_SNR(2024)
} else {
  all_checkers <- generate_checkers(2024)
}

print(paste("Getting files for path", path))

paths_to_combine <- list.files(path, full.names = TRUE)

run_ids <- as.integer(str_match(paths_to_combine, "_(\\d{1,3})\\.RDS")[, 2])

print(paste("Reading", length(paths_to_combine), "RDS files"))

all_rds <- list()

num_chunks <- 10

partition <- split(seq_along(paths_to_combine),
                   sample(1:num_chunks, length(paths_to_combine), replace=T))

as_csv <- tibble()
for(chunk in 1:num_chunks) {
  print(paste("Reading chunk", chunk))

  pb <- txtProgressBar(min = 0,
                       max = length(partition[[chunk]]),
                       initial = 0)

  chunk_rds <- list()
  t = 1
  for(path_index in partition[[chunk]]) {
    chunk_rds[[run_ids[path_index]]] <- readRDS(paths_to_combine[path_index])

    setTxtProgressBar(pb, t)
    t = t + 1
  }

  print(paste("Getting chunk", chunk, "results"))
  as_csv <- as_csv |>
    rbind(get_all_checker_results(all_checkers, chunk_rds, bcbc_extractor))

  print(paste("Garbage collecting", chunk))
  gc()
}

as_csv <- as_csv |> arrange(p, trial)

csv_path <- paste0("~/BCBC/simulated_experiments/rds_files/", pattern, ".csv")

print(paste("Saving csv to", csv_path))

write.csv(as_csv, csv_path, row.names=FALSE)
