#' Perform cross-validation with a hold-out test set. Used to find optimal gamma
#' hyperparameter.
#'
#' @param X data matrix
#' @param holdout_size percent of data to leave out
#' @param lambdas passed to `BCBC_missing`
#' @param gammas passed to `BCBC_missing`
#' @param k_samples passed to `BCBC_missing`
#' @param k_features passed to `BCBC_missing`
#' @param phis passed to `BCBC_missing`
#' @param recalculate_weights passed to `BCBC_missing`
#' @param approx_neighbors passed to `BCBC_missing`
#' @param hnsw_args passed to `BCBC_missing`
#' @param tols passed to `BCBC_missing`
#' @param tmax_hierarchy passed to `BCBC_missing`
#' @param return_fits passed to `BCBC_missing`
#' @param ... passed to `BCBC_missing`
#'
#' @return list with info
#'   1. `all_runs` all returned runs from `BCBC_missing`
#'   1. `gamma_cv_data` data frame of hyperparameters and cv metrics
#'   1. `all_params` grid of hyperparameters used for fitting
#'   1. `mask` of hold-out data
#'   1. `holdout_size`, `tmax_hierarchy` input parameters
#' @seealso [BCBC_missing()]
#' @seealso [cv.BCBC()]
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = FALSE)
#' gamma_cv <- cv.BCBC_holdout(
#'   checker$X,
#'   holdout_size = 0.2,
#'   lambdas = c(0),
#'   gammas = seq(50, 300, 50),
#'   k_samples = 4,
#'   k_features = 4,
#'   recalculate_weights = T
#' )
#' View(gamma_cv$gamma_cv_data)
cv.BCBC_holdout <- function(X,
                            holdout_size = 0.1,
                            lambdas = c(1),
                            gammas = c(1),
                            k_samples = c(2),
                            k_features = c(2),
                            phis = c(0.5),
                            recalculate_weights = c(FALSE),
                            approx_neighbors = c(FALSE),
                            hnsw_args = list(),
                            tols = c(1e-3),
                            tmax_hierarchy = c(10, 25, 50),
                            return_fits = FALSE,
                            ...) {
  all_params <-
    expand.grid(
      lambda = lambdas,
      k_samples = k_samples,
      k_features = k_features,
      gamma = gammas,
      phi = phis,
      recalculate_weights = recalculate_weights,
      approx_neighbors = approx_neighbors,
      tol = tols
    )

  cv_data <- data.frame(all_params)
  all_runs <- list()
  row_clusters <- list()
  col_clusters <- list()

  missing_X <- X
  mask <- which(matrix(runif(nrow(X) * ncol(X)), nrow = nrow(X)) < holdout_size, arr.ind = TRUE)
  heldout_vals <- X[mask]
  missing_X[mask] <- NA

  results <- lapply(1:nrow(all_params), function(param_set) {
    params <- all_params[param_set, ]

    # Use non-missing weights for CV
    wts <- fast_gkn_weights(
      t(X),
      k_row = params$k_features,
      k_col = params$k_samples,
      phi = params$phi,
      approximate = params$approx_neighbors,
      hnsw_args = hnsw_args
    )

    print(params)
    result <- do.call(BCBC_missing, c(list(X = missing_X,
                                           fusion_wts = wts,
                                           tmax_hierarchy = tmax_hierarchy,
                                           hnsw_args = hnsw_args,
                                           return_fits = return_fits),
                                      params, ...))

    residuals <- result$filled_vals[nrow(result$filled_vals), ] - heldout_vals
    cv_metrics <- list(rss_heldout = sum(residuals ^ 2))

    list(result = result, cv_metrics = cv_metrics)
  })

  for (param_set in 1:nrow(all_params)) {
    all_runs[[param_set]] = results[[param_set]]$result
    cv_metrics <- results[[param_set]]$cv_metrics
    cv_data[param_set, names(cv_metrics)] <- cv_metrics
  }

  return(list(
    all_runs = all_runs,
    gamma_cv_data = cv_data,
    all_params = all_params,
    mask = mask,
    holdout_size = holdout_size,
    tmax_hierarchy = tmax_hierarchy
  ))
}
