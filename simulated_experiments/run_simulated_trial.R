library(tidyverse)
library(BCBC)
library(future.apply)

args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 3)

JOB_ID <- args[1]
ARRAY_ID <- as.integer(args[2])
METHOD <- args[3]

plan(sequential, split=TRUE)

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

matched_algo <- match(toupper(METHOD), c("ADAPTIVE_BCBC", "BCBC", "BCBC_NO_SCALE",
                                         "BCBC_NO_SCALE_OR_ADAPT", "BCEL", "COBRA"))

if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0, 2 ^ seq(-7, 0, 1)),
    gammas = 2 ^ seq(1, 10, 2),
    k_rows = c(4),
    k_cols = c(10),
    phis = c(0.25),
    tols = c(1e-2),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_params <- result$cv_data |> slice(which.min(rss_heldout))
  print(best_params)

  result$bcbc <- BCBC(
    all_checkers[[ARRAY_ID]]$X,
    lambda = best_params$lambda,
    gamma = best_params$gamma,
    k_row = 4,
    k_col = 10,
    phi = 0.25,
    recalculate_weights = TRUE,
    scale_gamma = TRUE,
    progress = TRUE,
    tmax_outer = 500
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0, 2 ^ seq(-2, 6, 1)),
    gammas = 2 ^ seq(1, 10, 2),
    k_rows = c(4),
    k_cols = c(10),
    phis = c(0.25),
    tols = c(1e-2),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_params <- result$cv_data |> slice(which.min(rss_heldout))
  print(best_params)

  result$bcbc <- BCBC(
    all_checkers[[ARRAY_ID]]$X,
    lambda = best_params$lambda,
    gamma = best_params$gamma,
    k_row = 4,
    k_col = 6,
    phi = 0.25,
    recalculate_weights = FALSE,
    scale_gamma = TRUE,
    progress = TRUE,
    tmax_outer = 500
  )
} else if(matched_algo == 3) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0, 2 ^ seq(-7, -1, 1)),
    gammas = 2 ^ seq(3, 12, 1),
    k_rows = c(4),
    k_cols = c(6),
    phis = c(0.25),
    tols = c(10 ^ -2.5),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_params <- result$cv_data |> slice(which.min(rss_heldout))
  print(best_params)
  result$bcbc <- BCBC(
    all_checkers[[ARRAY_ID]]$X,
    lambda = best_params$lambda,
    gamma = best_params$gamma,
    k_row = 4,
    k_col = 6,
    phi = 0.25,
    tol = 10 ^ -4.5,
    recalculate_weights = TRUE,
    scale_gamma = FALSE,
    progress = TRUE,
    tmax_outer = 500
  )
} else if(matched_algo == 4) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0, 2 ^ seq(-7, -1, 1)),
    gammas = 2 ^ seq(3, 12, 1),
    k_rows = c(4),
    k_cols = c(6),
    phis = c(0.25),
    tols = c(10 ^ -2.5),
    tmax_outers = c(50),
    tmax_inners = c(50)
  )

  best_params <- result$cv_data |> slice(which.min(rss_heldout))
  print(best_params)

  result$bcbc <- BCBC(
    all_checkers[[ARRAY_ID]]$X,
    lambda = best_params$lambda,
    gamma = best_params$gamma,
    k_row = 4,
    k_col = 6,
    phi = 0.25,
    tol = 1e-4,
    recalculate_weights = FALSE,
    scale_gamma = FALSE,
    progress = TRUE,
    tmax_outer = 500
  )
} else if(matched_algo == 5) {
  library(BCEL)
  result <-
    bcel_stable(
      all_checkers[[ARRAY_ID]]$X,
      r = 5
    )
} else if(matched_algo == 6) {
  wts <- fast_gkn_weights(
    t(all_checkers[[ARRAY_ID]]$X),
    k_row = 2,
    k_col = 2,
    phi = 0.5,
    approx = 0
  )

  result <-
    cobra_validate(
      t(all_checkers[[ARRAY_ID]]$X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = 2 ^ seq(3, 16, 0.25),
      max_iter = 500
    )
} else {
  stop(c("Invalid args:", args))
}

saveRDS(result, paste0(
  "/work/sgr26/",
  toupper(METHOD),
  "_",
  JOB_ID,
  "_",
  ARRAY_ID,".RDS")
)

