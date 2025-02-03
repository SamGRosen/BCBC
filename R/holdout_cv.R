#' Title
#'
#' @param X
#' @param holdout_size
#' @param lambdas
#' @param gammas
#' @param k_samples
#' @param k_features
#' @param phis
#' @param recalculate_weights
#' @param tols
#' @param tmax_outers
#' @param tmax_inners
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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
