library(BCBC)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 3)

JOB_ID <- args[1]
METHOD <- args[2]
DATASET <- args[3]

matched_algo <- match(toupper(METHOD), c("BCBC_ADAPTIVE", "BCBC",
                                         "BCEL", "COBRA", "SPARSEBC",
                                         "SCBICLUST", "BCBC_ADAPT_APPROX"))

datasets <- c(
  "srbct.x.txt",
  "brain.x.txt",
  "colon.x.txt",
  "leukemia.x.txt",
  "lung100.csv",
  "lung500.csv",
  "lymphoma.x.txt",
  "prostate.x.txt",
  "WangBreast2005.txt",
  "BhattacherjeeLung2001.txt"
)
matched_data <- pmatch(
  DATASET,
  datasets
)

print(datasets)
print(DATASET)
print(matched_data)
X <- t(read.csv(paste0("./data/", datasets[matched_data]), header = FALSE, sep = " "))
X <- scale(X)

# BCBC params
k_samples = 10
k_features = 20
phi = 1
holdout_size = 0.25
tmax_hierarchy = c(2, 75, 300)
fit_tmax_hierarchy = c(4 * tmax_hierarchy[2], tmax_hierarchy[3])

p <- ncol(X)
n <- nrow(X)

gammas <- n * 2 ^ seq(-1, 6, 0.5)

lambdas <- seq(0, 200, 10) / p

debug <- TRUE


if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    X,
    holdout_size = holdout_size,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    recalculate_weights = c(TRUE),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = tmax_hierarchy,
    return_fits = debug
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$lambda_cv <- cv.BCBC(
    X,
    gammas = c(best_gamma$gamma),
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    percent_noise = seq(0.025, 0.50, 0.025),
    recalculate_weights = c(TRUE),
    tols = c(1e-4),
    tmax_outer = fit_tmax_hierarchy[1],
    tmax_cobra = fit_tmax_hierarchy[2]
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    X,
    holdout_size = holdout_size,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    recalculate_weights = c(FALSE),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = c(2, 150, 250),
    return_fits = debug
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$lambda_cv <- cv.BCBC(
    X,
    gammas = c(best_gamma$gamma),
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    percent_noise = seq(0.025, 0.50, 0.025),
    recalculate_weights = c(FALSE),
    tols = c(1e-4),
    tmax_outer = fit_tmax_hierarchy[1],
    tmax_cobra = fit_tmax_hierarchy[2]
  )
} else if (matched_algo == 3) {
  library(BCEL)
  # TODO Need to give this true number of clusters
  result <- bcel_stable(X, r = 5)
} else if(matched_algo == 4) {
  wts <- fast_gkn_weights(
    t(X),
    k_row = k_features,
    k_col = k_samples,
    phi = phi,
    approx = 0
  )

  cobra_gammas <- 2 ^ seq(3, 16, 0.25)
  gamma_cv <-
    cobra_validate(
      t(X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = cobra_gammas,
      max_iter = 500
    )

  best_gamma <- cobra_gammas[which.min(gamma_cv$validation_error)]

  result <-
    cobra(
      t(X),
      wts$E_row,
      wts$E_col,
      wts$w_row,
      wts$w_col,
      gamma = best_gamma,
      max_iter = 500
    )
} else if(matched_algo == 5) {
  library(sparseBC)
  sparsebc_gammas <- 2 ^ seq(3, 16, 0.25)
  tuning_start <- Sys.time()
  kr <-
    sparseBC.choosekr(
      X,
      k = 2:5,
      r = seq(5, 50, 5),
      lambda = 0,
      percent = 0.125,
      TRUE
    )

  print(kr)

  if(is.null(kr$estimated_kr)) {
    best_k <- kr$bestK
    best_r <- kr$bestR
  } else {
    best_k <- kr$estimated_kr[1, 1]
    best_r <- kr$estimated_kr[1, 2]
  }

  lambda <- sparseBC.BIC(X,
                         best_k,
                         best_r,
                         sparsebc_gammas)
  print(paste(
    "Tuning time:",
    as.numeric(Sys.time() - tuning_start, units = "secs")
  ))

  fit_start <- Sys.time()

  result <- sparseBC(X,
                     best_k,
                     best_r,
                     lambda$lambda)
  print(paste("Fit time:", as.numeric(Sys.time() - fit_start, units = "secs")))
  print(warnings())
} else if(matched_algo == 6) {
  library(SCBiclust)
  result <- PermBiclust.sigclust_stop(X,
                                      sc = FALSE,
                                      silent = TRUE,
                                      maxnum.bicluster = 5)
} else if(matched_algo == 7) {
  result <- cv.BCBC_holdout(
    X,
    holdout_size = holdout_size,
    lambdas = c(0),
    gammas = gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -4),
    tmax_hierarchy = tmax_hierarchy,
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
    X,
    gammas = c(best_gamma$gamma),
    lambdas = lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(TRUE),
    tols = c(1e-4),
    tmax_outer = fit_tmax_hierarchy[1],
    tmax_cobra = fit_tmax_hierarchy[2],
    approx_neighbors = c(TRUE),
    percent_noise = seq(0.025, 0.50, 0.025),
    hnsw_args = list(
      ef = 100,
      M = 32,
      n_threads = 4
    )
  )
} else {
  stop(c("Invalid args:", args))
}


output_dir <- "/cwork/sgr26/"

saveRDS(result, paste0(
  output_dir,
  toupper(METHOD),
  "_",
  DATASET,
  "_",
  JOB_ID,
  ".RDS")
)
