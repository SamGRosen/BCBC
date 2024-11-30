library(BCBC)
library(dplyr)
library(future.apply)

args <- commandArgs(trailingOnly = TRUE)

print(args)

stopifnot(length(args) == 3)

JOB_ID <- args[1]
METHOD <- args[2]
DATASET <- args[3]

plan(sequential, split=TRUE)

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

X <- t(read.csv(paste0("./data/", datasets[matched_data]), header = FALSE, sep = " "))
X <- scale(X)


k_samples = 4
k_features = 10
phi = 0.25

p <- ncol(X)

bcbc_gammas <- p * 2 ^ seq(-2, 1, 0.5)  # change back to more, saving time in cv
bcbc_gammas <- p * 2 ^ seq(0, 1, 0.25)  # change back to more, saving time in cv

bcbc_lambdas <- c(0, 2 ^ seq(-9, -1, 1))

if(matched_algo == 1) {
  result <- cv.BCBC_holdout(
    X,
    holdout_size = 0.2,
    lambdas = c(0),
    gammas = bcbc_gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    recalculate_weights = c(TRUE),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50)
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$bcbc_cv <- cv.BCBC(
    X,
    gammas = c(best_gamma$gamma),
    lambdas = bcbc_lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(TRUE),
    tols = c(1e-4),
    tmax_outers = c(500)
  )
} else if(matched_algo == 2) {
  result <- cv.BCBC_holdout(
    X,
    holdout_size = 0.15,
    lambdas = c(0),
    gammas = bcbc_gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(15, 150, 50)
  )

  best_gamma <- result$gamma_cv_data |> slice(which.min(rss_heldout))

  print(best_gamma)

  result$bcbc_cv <- cv.BCBC(
    X,
    gammas = c(best_gamma$gamma),
    lambdas = bcbc_lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(FALSE),
    tols = c(1e-4),
    tmax_outers = c(500)
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
    holdout_size = 0.10,
    lambdas = c(0),
    gammas = bcbc_gammas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    tols = c(10 ^ -2.5),
    tmax_hierarchy = c(10, 50, 50),
    recalculate_weights = c(TRUE),
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
    X,
    gammas = c(best_gamma$gamma),
    lambdas = bcbc_lambdas,
    k_samples = c(k_samples),
    k_features = c(k_features),
    phis = c(phi),
    recalculate_weights = c(TRUE),
    tols = c(1e-4),
    tmax_outers = c(500),
    approx_neighbors = c(TRUE),
    hnsw_args = list(
      ef = 100,
      M = 32,
      n_threads = 12
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
