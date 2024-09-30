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
#' @import future.apply
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
                            tols = c(1e-3),
                            tmax_outers = c(50),
                            tmax_inners = c(50),
                            ...) {
  all_params <-
    expand.grid(
      lambda = lambdas,
      k_samples = k_samples,
      k_features = k_features,
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
  mask <- which(matrix(runif(nrow(X) * ncol(X)), nrow = nrow(X)) < holdout_size, arr.ind = TRUE)
  heldout_vals <- X[mask]
  missing_X[mask] <- NA

  future_results <- lapply(1:nrow(all_params), function(param_set) {
    params <- all_params[param_set, ]
    print(params)
    result <- do.call(BCBC_missing, c(list(X = missing_X), params, ...))

    residuals <- result$filled_vals[nrow(result$filled_vals), ] - heldout_vals
    unique_weights <- result$w_path[nrow(result$w_path), ]

    weights <- result$w_path[nrow(result$w_path), mask[, 2]]
    nonzeros <- which(weights > 0)
    cv_metrics <- list(rss_heldout = sum(residuals ^ 2))

    list(result = result, cv_metrics = cv_metrics)
  })

  for (param_set in 1:nrow(all_params)) {
    all_runs[[param_set]] = future_results[[param_set]]$result
    cv_metrics <- future_results[[param_set]]$cv_metrics
    cv_data[param_set, names(cv_metrics)] <- cv_metrics
  }

  return(list(
    all_runs = all_runs,
    gamma_cv_data = cv_data,
    all_params = all_params,
    mask = mask
  ))
}


#' An alternative method to select a value for lambda by checking if the clusters
#' from BCBC match the clusters from COBRA run on the weighted matrix.
#'
#' @param X
#' @param lambdas
#' @param gamma
#' @param k_samples
#' @param k_features
#' @param phi
#' @param tol
#' @param tmax_outer
#' @param tmax_inner
#' @param reweight reweight the columns before passing into COBRA
#' @param percent_noise
#' @param progress
#' @param recalculate_weights
#'
#' @return
#' @export
#'
#' @examples
cv.BCBC_holdout_stability <- function(X,
                                      gamma = 1,
                                      lambdas = c(0),
                                      k_samples = 2,
                                      k_features = 2,
                                      phi = 0.5,
                                      recalculate_weights = TRUE,
                                      tol = 1e-4,
                                      tmax_outer = 50,
                                      tmax_inner = 50,
                                      reweight = TRUE,
                                      percent_noise = 0.25,
                                      method = "adjusted.rand",
                                      progress = TRUE) {
  lambda_fits <- list()
  cv_data <- data.frame(
    lambda = numeric(),
    gamma = numeric(),
    row_similarity = numeric(),
    col_similarity = numeric(),
    num_cobra_row_clusters = numeric(),
    num_cobra_col_clusters = numeric(),
    num_bcbc_row_clusters = numeric(),
    num_bcbc_col_clusters = numeric()
  )
  i <- 1
  for (lambda in lambdas) {
    bcbc_fit <- BCBC(
      X,
      lambda = lambda,
      gamma = gamma,
      k_samples = k_samples,
      k_features = k_features,
      phi = phi,
      tmax_cobra = tmax_inner,
      tmax_outer = tmax_outer,
      recalculate_weights = recalculate_weights,
      progress = progress
    )
    lambda_fits[[i]] <- bcbc_fit

    non_zero <- which(bcbc_fit$w > 0)
    non_zero_w <- bcbc_fit$w[non_zero]

    bcbc_biclusters <- get_weighted_biclusters(bcbc_fit$U,
                                               weights = bcbc_fit$w,
                                               lambda = lambda,
                                               percent_noise = percent_noise)
    if (reweight) {
      to_cobra <- sweep(X[, non_zero], 2, non_zero_w ^ 2 + lambda * non_zero_w, "*")
    } else {
      to_cobra <- X[, non_zero]
    }

    cobra_fit <- cobra_knn(
      to_cobra,
      k_samples = k_row,
      k_features = k_col,
      phi = phi,
      tol = tol,
      gamma = best_gamma
    )

    cobra_biclusters <- get_weighted_biclusters(cobra_fit,
                                                weights = 1,
                                                percent_noise = percent_noise)

    row_similarity <- compare_clusters(cobra_biclusters$row_clusters,
                                       bcbc_biclusters$row_clusters,
                                       method = method)

    col_similarity <- compare_clusters(cobra_biclusters$col_clusters,
                                       bcbc_biclusters$col_clusters[non_zero_w],
                                       method = method)

    cv_data <- cv_data |>
      add_row(
        row_similarity = row_similarity,
        col_similarity = col_similarity,
        num_cobra_row_clusters = length(unique(cobra_biclusters$row_clusters)),
        num_cobra_col_clusters = length(unique(cobra_biclusters$col_clusters)),
        num_bcbc_row_clusters = length(unique(bcbc_biclusters$row_clusters)),
        num_bcbc_col_clusters = length(unique(bcbc_biclusters$col_clusters[non_zero_w])),
        lambda = lambda,
        gamma = gamma
      )

    i <- i + 1
  }

  return(list(
    lambda_fits = lambda_fits,
    cv_data = cv_data,
    mask = mask
  ))
}


compare_clusters <- function(clusters1, clusters2, method = "adjusted.rand") {
  if (length(unique(clusters1)) == length(clusters1) ||
      length(unique(clusters2)) == length(clusters2)) {
    return(-1)
  }
  igraph::compare(clusters1, clusters2, method = method)
}
