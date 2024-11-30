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

all_checkers <- generate_checkers_SNR2(2024)

if(ARRAY_ID > length(all_checkers)) {
  stop(paste("Array ID is larger than problem ID", ARRAY_ID, length(all_checkers)))
}

matched_algo <- match(toupper(METHOD), c("BCBC_NO_SCALE", "BCBC_NO_SCALE_OR_ADAPT",
                                         "BCBC_ADAPT_APPROX", "BCBC_ORACLE",
                                         "BCEL", "COBRA"))

k_samples <- 5
k_features <- 5
phi <- 1

p <- ncol(all_checkers[[ARRAY_ID]]$X)
if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.30,
    lambdas = c(0),
    gammas = p * 2 ^ seq(-2, 1, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50),
    recalculate_weights = c(TRUE)
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
    tols = c(1e-5),
    tmax_outers = c(300),
    tmax_cobras = c(100)
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.30,
    lambdas = c(0),
    # gammas = c(p * sqrt(2)),
    gammas = p * 2 ^ seq(-2, 1, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50),
    recalculate_weights = c(FALSE)
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$bcbc_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gammas = c(best_gamma$gamma),
    # lambdas = c(0),
    lambdas = c(0, 2 ^ seq(-7, -1, 1)),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(FALSE),
    greedy_terminate = FALSE,
    tols = c(1e-5),
    tmax_outers = c(300), tmax_cobras = c(100)
  )
} else if(matched_algo == 3) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.30,
    lambdas = c(0),
    gammas = p * 2 ^ seq(-2, 1, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50),
    recalculate_weights = c(FALSE)
    # approx_neighbors = c(TRUE),
    # hnsw_args = list(
    #   ef = 50,
    #   M = 32,
    #   n_threads = 4
    # )
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
    tols = c(1e-5),
    approx_neighbors = c(TRUE),
    hnsw_args = list(
      ef = 50,
      M = 16,
      n_threads = 4
    ),
    tmax_outers = c(300),
    tmax_cobras = c(100)
  )
} else if(matched_algo == 4) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.15,
    lambdas = c(0),
    gammas = p * 2 ^ seq(-3, 2, 0.5),
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50),
    recalculate_weights = c(FALSE)
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
  )
} else if (matched_algo == 5) {
  library(BCEL)
  result <- bcel_stable(all_checkers[[ARRAY_ID]]$X, r = 5)
} else if(matched_algo == 6) {
  wts <- fast_gkn_weights(
    t(all_checkers[[ARRAY_ID]]$X),
    k_row = k_features,
    k_col = k_samples,
    phi = phi,
    approx = 0
  )

  gammas <- 2 ^ seq(3, 15.5, 0.25)
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
  print(paste("Creating dir", output_dir))
  dir.create(output_dir)
} else {
  print(paste("Saving to existing", output_dir))
}



save_success <- tryCatch(
  expr = {
    saveRDS(result,
            paste0(
              output_dir,
              toupper(METHOD),
              "_",
              JOB_ID,
              "_",
              ARRAY_ID,
              ".RDS"
            ))
    TRUE
  },
  error = function(e) {
    print(e)
    FALSE
  }
)

if (!save_success) {
  q(save = "no",
    status = 11,
    runLast = FALSE)
}
