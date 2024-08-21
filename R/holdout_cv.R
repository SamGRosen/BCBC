#' Title
#'
#' @param X
#' @param holdout_size
#' @param lambdas
#' @param gammas
#' @param k_rows
#' @param k_cols
#' @param phis
#' @param recalculate_weights
#' @param tols
#' @param tmax_outers
#' @param tmax_inners
#' @param ...
#'
#' @return
#' @import future.apply
#' @export
#'
#' @examples
cv.BCBC_holdout <- function(X,
                            holdout_size = 0.1,
                            lambdas = c(1),
                            gammas = c(1),
                            k_rows = c(2),
                            k_cols = c(2),
                            phis = c(0.5),
                            recalculate_weights = c(FALSE),
                            tols = c(1e-3),
                            tmax_outers = c(50),
                            tmax_inners = c(50),
                            retrain = TRUE,
                            ...) {


  all_params <-
    expand.grid(
      lambda = lambdas,
      k_row = k_rows,
      k_col = k_cols,
      gamma = gammas,
      tmax_outer = tmax_outers,
      tmax_inner = tmax_inners,
      phi = phis,
      recalculate_weights = recalculate_weights,
      tol = tols
    )

  cv_data <- data.frame(all_params)
  all_runs <- list()
  row_clusters <- list()
  col_clusters <- list()

  missing_X <- X
  mask <- which(matrix(runif(nrow(X) * ncol(X)), nrow=nrow(X)) < holdout_size, arr.ind = TRUE)
  heldout_vals <- X[mask]
  missing_X[mask] <- NA

  future_results <- lapply(
    1:nrow(all_params),
    # future.seed = TRUE,
    function(param_set) {
      params <- all_params[param_set,]
      print(params)
      result <- do.call(BCBC_missing,
                        c(list(X=missing_X), params, ...))

      residuals <- result$filled_vals[nrow(result$filled_vals), ] - heldout_vals
      unique_weights <- result$w_path[nrow(result$w_path), ]

      weights <- result$w_path[nrow(result$w_path), mask[, 2]]
      sparsity_kurtosis <- sum(unique_weights ^ 2) ^ 2 / sum(unique_weights ^ 4)
      sparsity_kurtosis2 <- sum(weights ^ 2) ^ 2 / sum(weights ^ 4)

      nonzeros <- which(weights > 0)
      cv_metrics <- list(
        kurtosis = sparsity_kurtosis,
        total_kurtosis = sparsity_kurtosis2,
        rss_heldout = sum(residuals ^ 2),
        rss_heldout_abs = sum(abs(residuals)),
        rss_heldout_weighted = sum(residuals ^ 2 * weights),
        rss_heldout_weighted_sq = sum(residuals ^ 2 * weights ^ 2),
        rss_heldout_weighted_abs = sum(abs(residuals) * weights),
        rss_heldout_mean = mean(residuals[nonzeros]^2),
        rss_heldout_weighted_mean = mean(residuals[nonzeros]^2 * weights[nonzeros]),
        rss_heldout_weighted_mean_abs = mean(abs(residuals[nonzeros]) * weights[nonzeros]),
        rss_heldout_kurtosis = sum(residuals ^ 2 * weights) / sparsity_kurtosis,
        rss_heldout_kurtosis2 = sum(residuals ^ 2) / sparsity_kurtosis,
        rss_heldout_kurtosis3 = sum(residuals ^ 2 * weights) / sparsity_kurtosis2,
        rss_heldout_kurtosis4 = sum(residuals ^ 2) / sparsity_kurtosis2
      )

      list(
        result = result,
        cv_metrics = cv_metrics
      )
    })

  for(param_set in 1:nrow(all_params)) {
    all_runs[[param_set]] = future_results[[param_set]]$result
    cv_metrics <- future_results[[param_set]]$cv_metrics
    cv_data[param_set, names(cv_metrics)] <- cv_metrics
  }

  return(list(all_runs = all_runs,
              cv_data = cv_data,
              all_params = all_params,
              mask = mask))
}
