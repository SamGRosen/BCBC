library(dplyr)
library(BCBC)
source("~/BCBC/simulated_experiments/simulation_helpers.R")

options(warn = 1)


args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 3)

JOB_ID <- args[1]
ARRAY_ID <- as.integer(args[2])
METHOD <- args[3]

all_checkers <- generate_checkers(2024)

if(is.null(all_checkers[[ARRAY_ID]])) {
  stop(paste("Array ID is invalid for checkers", ARRAY_ID, length(all_checkers)))
}

matched_algo <- match(toupper(METHOD), c("BCBC_NO_SCALE", "BCBC_NO_SCALE_OR_ADAPT",
                                         "BCBC_ADAPT_APPROX", "BCEL", "COBRA",
                                         "COBRA_PREPROCESS"))

n <- nrow(all_checkers[[ARRAY_ID]]$X)
p <- ncol(all_checkers[[ARRAY_ID]]$X)
k_samples <- 5
k_features <- 5
phi <- 1
lambdas <- seq(0, 160, 4) / p
gammas = n * 1.5 ^ seq(-4, 4, 1)

holdout_out_tmax_hierarchy = c(2, 200, 200)

fit_tmax_hierarchy = holdout_out_tmax_hierarchy[2:3]

debug = TRUE

if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.25,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = holdout_out_tmax_hierarchy,
    recalculate_weights = c(TRUE),
    return_fits = debug
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$lambda_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gamma = c(best_gamma$gamma),
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phi = c(phi),
    percent_noise = seq(0.025, 0.25, 0.025),
    recalculate_weights = TRUE,
    tols = c(1e-4),
    tmax_outer = fit_tmax_hierarchy[1],
    tmax_cobra = fit_tmax_hierarchy[2]
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.25,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = c(2, 50, 200), # holdout_out_tmax_hierarchy,
    recalculate_weights = c(FALSE),
    return_fits = debug
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$lambda_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gamma = best_gamma$gamma,
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phi = c(phi),
    recalculate_weights = c(FALSE),
    percent_noise = seq(0.025, 0.25, 0.025),
    tols = c(1e-4),
    tmax_outer = 300, # fit_tmax_hierarchy[1],
    tmax_cobra = 200 # fit_tmax_hierarchy[2],
  )
} else if(matched_algo == 3) {
  result <- cv.BCBC_holdout(
    all_checkers[[ARRAY_ID]]$X,
    holdout_size = 0.25,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = holdout_out_tmax_hierarchy,
    recalculate_weights = c(TRUE),
    approx_neighbors = c(TRUE),
    hnsw_args = list(
      ef = 50,
      M = 32,
      n_threads = 4
    ),
    return_fits = debug
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$lambda_cv <- cv.BCBC(
    all_checkers[[ARRAY_ID]]$X,
    gammas = c(best_gamma$gamma),
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(TRUE),
    tols = c(1e-4),
    approx_neighbors = c(TRUE),
    percent_noise = seq(0.025, 0.25, 0.025),
    hnsw_args = list(
      ef = 50,
      M = 16,
      n_threads = 4
    ),
    tmax_outer = fit_tmax_hierarchy[1],
    tmax_cobra = fit_tmax_hierarchy[2]
  )
} else if (matched_algo == 4) {
  library(BCEL)
  result <- bcel_stable(all_checkers[[ARRAY_ID]]$X, r = 5)
} else if(matched_algo == 5) {
  wts <- fast_gkn_weights(
    t(all_checkers[[ARRAY_ID]]$X),
    k_row = k_features,
    k_col = k_samples,
    phi = phi,
    approx = 0
  )

  gammas <- 2 ^ seq(-3, 15, 0.25)
  gamma_cv <-
    cobra_validate(
      t(all_checkers[[ARRAY_ID]]$X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = gammas,
      max_iter = 2000,
      fraction = 0.25
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
      max_iter = 2000
    )
} else if(matched_algo == 6) {
  gammas <- 2 ^ seq(-3, 15, 0.25)

  top_k_indices <- function(vec, k) {
    order(vec, decreasing = TRUE)[1:k]
  }

  scales <- attr(all_checkers[[ARRAY_ID]]$X, "scaled:scale")
  top_200 <- top_k_indices(scales, 200)

  X_prime <- all_checkers[[ARRAY_ID]]$X[, top_200]
  wts <- fast_gkn_weights(
    t(X_prime),
    k_row = k_features,
    k_col = k_samples,
    phi = phi,
    approx = 0
  )

  gamma_cv <-
    cobra_validate(
      t(X_prime),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = gammas,
      max_iter = 2000,
      fraction = 0.25
    )

  best_gamma <- gammas[which.min(gamma_cv$validation_error)]

  result <-
    cobra(
      t(X_prime),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = best_gamma,
      max_iter = 2000
    )
  result$top_200 <- top_200
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
