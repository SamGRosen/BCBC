library(tidyverse)
library(BCBC)
library(future.apply)
source("~/BCBC/simulated_experiments/simulation_helpers.R")

args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 3)

JOB_ID <- args[1]
ARRAY_ID <- as.integer(args[2])
METHOD <- args[3]

plan(sequential, split=TRUE)

all_checkers <- generate_checkers(2024)

if(ARRAY_ID > length(all_checkers)) {
  stop(paste("Array ID is larger than problem ID", ARRAY_ID, length(all_checkers)))
}

matched_algo <- match(toupper(METHOD), c("BCBC_NO_SCALE", "BCBC_NO_SCALE_OR_ADAPT",
                                         "BCEL", "COBRA"))

k_samples <- 5
k_features <- 5
phi <- 1

p <- ncol(all_checkers[[ARRAY_ID]]$X)
if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0),
    gammas = p * 2 ^ seq(-3, 2, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$bcbc_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gammas = c(best_gamma$gamma),
    lambdas = c(0, 2 ^ seq(-7, -1, 1)),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(TRUE),
    tols = c(1e-5)
    # num_row_clusters = 5,
    # num_col_clusters = 6
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0),
    gammas = p * 2 ^ seq(-3, 2, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$bcbc_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gammas = c(best_gamma$gamma),
    lambdas = c(0, 2 ^ seq(-7, -1, 1)),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(FALSE),
    tols = c(1e-5)
    # num_row_clusters = 5,
    # num_col_clusters = 6
  )
} else if (matched_algo == 3) {
  library(BCEL)
  result <- bcel_stable(all_checkers[[ARRAY_ID]]$X, r = 5)
} else if(matched_algo == 4) {
  wts <- fast_gkn_weights(
    t(all_checkers[[ARRAY_ID]]$X),
    k_row = k_features,
    k_col = k_samples,
    phi = phi,
    approx = 0
  )

  gammas <- 2 ^ seq(3, 16, 0.25)
  gamma_cv <-
    cobra_validate(
      t(all_checkers[[ARRAY_ID]]$X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = gammas,
      max_iter = 500
    )

  best_gamma <- gammas[which.min(gamma_cv$validation_error)]

  result <-
    cobra(
      t(all_checkers[[ARRAY_ID]]$X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = best_gamma,
      max_iter = 500
    )
} else {
  stop(c("Invalid args:", args))
}


output_dir <- paste0(
  "/cwork/sgr26/",
  toupper(METHOD),
  "_",
  JOB_ID,
  "/"
)

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

saveRDS(result, paste0(
  output_dir,
  toupper(METHOD),
  "_",
  JOB_ID,
  "_",
  ARRAY_ID,
  ".RDS")
)

