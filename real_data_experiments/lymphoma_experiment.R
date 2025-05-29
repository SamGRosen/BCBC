library(BCBC)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2)

JOB_ID <- args[1]
METHOD <- args[2]

matched_algo <- match(toupper(METHOD), c("BCBC_ADAPTIVE", "BCBC",
                                         "BCBC_ADAPT_APPROX"))

X <- t(read.csv("./data/lymphoma_labelled.txt", header = TRUE, sep = ","))
X <- scale(X)

# BCBC params
k_samples = 10
k_features = 25
phi = 1
tmax_hierarchy = c(5, 75, 300)
fit_tmax_hierarchy = c(6 * tmax_hierarchy[2], tmax_hierarchy[3])

p <- ncol(X)
n <- nrow(X)

global_gamma <- 750

lambdas <- seq(0, 100, 5) / p

debug <- TRUE

if(matched_algo == 1) {
  result <- BCBC_missing(
    X,
    lambda = 0,
    k_samples = k_samples,
    k_features = k_features,
    gamma = global_gamma,
    phi = phi,
    tmax_hierarchy = tmax_hierarchy,
    tol = 1e-10,
    recalculate_weights = TRUE,
    return_fits = debug
  )

  filled_X = X
  filled_X[result$invalid_indices] = result$filled_vals[tmax_hierarchy[1],]

  result$filled_X = filled_X

  result$lambda_cv <- cv.BCBC(
    filled_X,
    gammas = global_gamma,
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
  result <- BCBC_missing(
    X,
    lambda = 0,
    k_samples = k_samples,
    k_features = k_features,
    gamma = global_gamma,
    phi = phi,
    tmax_hierarchy = tmax_hierarchy,
    tol = 1e-10,
    recalculate_weights = FALSE,
    return_fits = debug
  )

  filled_X = X
  filled_X[result$invalid_indices] = result$filled_vals[tmax_hierarchy[1],]

  result$lambda_cv <- cv.BCBC(
    filled_X,
    gammas = global_gamma,
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
  result <- BCBC_missing(
    X,
    lambda = 0,
    k_samples = k_samples,
    k_features = k_features,
    gamma = global_gamma,
    phi = phi,
    tmax_hierarchy = tmax_hierarchy,
    tol = 1e-10,
    recalculate_weights = TRUE,
    return_fits = debug,
    approx_neighbors = TRUE,
    hnsw_args = list(
      ef = 50,
      M = 32,
      n_threads = 4
    )
  )

  filled_X = X
  filled_X[result$invalid_indices] = result$filled_vals[tmax_hierarchy[1],]

  result$lambda_cv <- cv.BCBC(
    filled_X,
    gammas = global_gamma,
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
  "lymphoma_labelled",
  "_",
  JOB_ID,
  ".RDS")
)
