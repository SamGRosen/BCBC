#' Get metrics to measure a single BCBC fit that are useful for cross validation
#'
#' @param X data matrix
#' @param bcbc_run return value from `BCBC`
#' @param percent_noise Percent of standard deviation of pairwise distances to use for thresholding
#' @param num_row_clusters If set, tries to calculate a threshold to result in this many row clusters
#' @param num_col_clusters If set, tries to calculate a threshold to result in this many column clusters
#'
#' @return list with
#'   1. `cv_metrics` data frame with cv info such as RSS and eBIC
#'   2. `row_clusters` calculated row clusters
#'   3. `col_clusters` calculated column clusters
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = FALSE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.05, gamma = 200)
#' get_cv_metrics(checker$X, bcbc_fit, percent_noise = 0.25)
get_cv_metrics <- function(X,
                           bcbc_run,
                           percent_noise = 0.25,
                           num_row_clusters = NA,
                           num_col_clusters = NA) {
  n <- nrow(X)
  p <- ncol(X)

  biclusters <- bicluster_assignments(
    bcbc_run$U,
    weights = bcbc_run$w,
    lambda = bcbc_run$lambda,
    num_row_clusters = num_row_clusters,
    num_col_clusters = num_col_clusters,
    percent_noise = percent_noise
  )

  fitted_biclust = bicluster_centers(X, biclusters$row_clusters, biclusters$col_clusters)

  U_rss <- sum((X - bcbc_run$U)^2)
  bicluster_rss <- sum((X - fitted_biclust$fit)^2)
  w2_rss <- sum(bcbc_run$w ^ 2 * colSums((X - bcbc_run$U)^2))
  w0 <- sum(bcbc_run$w > 0)
  row_clusters <- biclusters$row_clusters
  col_clusters <- biclusters$col_clusters

  num_row_clusters <- length(unique(row_clusters))
  num_col_clusters <- length(unique(col_clusters)) - (w0 != p)
  num_bi_clusters <- length(unique(c(fitted_biclust$assignments)))

  list(
    cv_metrics = list(
      num_row_clusters = num_row_clusters,
      num_col_clusters = num_col_clusters,
      U_rss = U_rss,
      bicluster_rss = bicluster_rss,
      w2_rss = w2_rss,
      BIC = n * p * log(bicluster_rss / (n * p)) + log(n * p) * num_bi_clusters,
      eBIC = n * p * log(bicluster_rss / (n * p)) + 2 * log(n * p) * num_bi_clusters,
      eBIC2 = n * p * log(bicluster_rss / (n * p)) + 3 * log(n * p) * num_bi_clusters,
      w0 = w0,
      percent_noise = percent_noise
    ),
    row_clusters = row_clusters,
    col_clusters = col_clusters
  )
}


#' Utility method to get cv metrics for many BCBC fits. Use `cv.BCBC` instead.
#'
#' @param X data matrix
#' @param bcbc_runs list of bcbc runs to use
#' @param all_params grid of hyperparameters
#' @param percent_noise vector of potential percent_noises to build biclusters
#' @inheritParams get_cv_metrics
#' @seealso [cv.BCBC()]
#'
#' @return list with
#'   1. `all_runs` used for calculation
#'   1. `cv_data` all metrics used for cv
#'   1. `row_clusters` all row clusters
#'   1. `col_clusters` all column clusters
#'   1. `all_params` hyperparameter grid
#' @export
get_all_cv_metrics <- function(X,
                               bcbc_runs,
                               all_params = data.frame(),
                               percent_noise = c(0.25),
                               num_row_clusters = NA,
                               num_col_clusters = NA) {
  cv_data <- data.frame()
  row_clusters <- list()
  col_clusters <- list()

  index <- 1
  for (param_index in 1:length(bcbc_runs)) {
    for (percent in percent_noise) {
      cv_info <- get_cv_metrics(
        X,
        bcbc_runs[[param_index]],
        percent_noise = percent,
        num_row_clusters = num_row_clusters,
        num_col_clusters = num_col_clusters
      )
      cv_data[index, "param_index"] <- param_index
      cv_data[index, names(all_params)] <- all_params[param_index,]
      cv_data[index, names(cv_info$cv_metrics)] <- cv_info$cv_metrics
      row_clusters[[index]] <- cv_info$row_clusters
      col_clusters[[index]] <- cv_info$col_clusters

      index <- index + 1
    }
  }

  cv_data <- cv_data |>
    mutate(index = row_number())

  return(
    list(
      all_runs = bcbc_runs,
      cv_data = cv_data,
      row_clusters = row_clusters,
      col_clusters = col_clusters,
      all_params = all_params
    )
  )
}

#' Perform many runs of `BCBC` for cross-validation purposes.
#'
#' @param X data matrix
#' @param lambdas inputs to `BCBC`
#' @param k_samples inputs to `BCBC`
#' @param k_features inputs to `BCBC`
#' @param gammas inputs to `BCBC`
#' @param tmaxs inputs to `BCBC`
#' @param tmax_outers inputs to `BCBC`
#' @param tmax_cobras inputs to `BCBC`
#' @param phis inputs to `BCBC`
#' @param recalculate_weights inputs to `BCBC`
#' @param tols inputs to `BCBC`
#' @param percent_noise input to `get_cv_metrics`
#' @param progress show progress
#' @param ... passed to BCBC
#'
#' @inherit get_all_cv_metrics return
#' @seealso [BCBC()]
#' @seealso [cv.BCBC_holdout()]
#' @seealso [get_cv_metrics()]
#' @seealso [get_all_cv_metrics()]
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit_info <- cv.BCBC(
#'   checker$X,
#'   lambdas = seq(0, 0.2, 0.05),
#'   percent_noise = seq(0.05, 0.25, 0.05),
#'   gammas = 100
#' )
#'
#' best_run_info <- bcbc_fit_info$cv_data |>
#'   filter(is.finite(eBIC2)) |>
#'   slice(which.min(eBIC2))
#' best_param_index <- best_run_info$param_index
#' best_cluster_index <- best_run_info$index
#' row_groups <- bcbc_fit_info$row_clusters[[best_cluster_index]]
#' col_groups <- bcbc_fit_info$col_clusters[[best_cluster_index]]
#' best_run <- bcbc_fit_info$all_runs[[best_param_index]]
#'
#' fit_as_df <- bcbc_result_to_df(best_run, labels = row_groups, filter_weight = 0)
#' plot_fit(fit_as_df, bin_scale = 10)
cv.BCBC <- function(X,
                    lambdas = c(1),
                    k_samples = c(2),
                    k_features = c(2),
                    gammas = c(1),
                    tmaxs = c(NA),
                    tmax_outers = c(100),
                    tmax_cobras = c(100),
                    phis = c(0.5),
                    recalculate_weights = c(TRUE),
                    tols = c(1e-6),
                    percent_noise = c(0.25),
                    approx_neighbors = c(FALSE),
                    hnsw_args = list(),
                    num_row_clusters = NA,
                    num_col_clusters = NA,
                    ...) {

  all_params <-
    expand.grid(
      lambda = lambdas,
      k_samples = k_samples,
      k_features = k_features,
      gamma = gammas,
      tmax = tmaxs,
      tmax_outer = tmax_outers,
      tmax_cobra = tmax_cobras,
      phi = phis,
      recalculate_weights = recalculate_weights,
      approx_neighbors = approx_neighbors,
      tol = tols
    )

  cv_data <- data.frame(all_params) |>
    mutate(index=row_number())

  all_runs <- list()
  row_clusters <- list()
  col_clusters <- list()

  for(param_index in 1:nrow(all_params)) {
    params <- all_params[param_index,]
    all_runs[[param_index]] = do.call(BCBC,
                                      c(list(X=X, hnsw_args=hnsw_args),
                                        params,
                                        ...))
  }

  get_all_cv_metrics(
    X,
    all_runs,
    all_params,
    percent_noise = percent_noise,
    num_row_clusters = num_row_clusters,
    num_col_clusters = num_col_clusters
  )
}
