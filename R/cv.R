#' Get metrics to measure a single BCBC fit that are useful for cross validation
#'
#' @param X
#' @param bcbc_run
#' @param weighted_clusters
#' @param percent_noise Percent of standard deviation of pairwise distances to use for thresholding
#' @param num_row_clusters If set, tries to calculate a threshold to result in this many row clusters
#' @param num_col_clusters If set, tries to calculate a threshold to result in this many column clusters
#'
#' @return
#' @export
#'
#' @examples
get_cv_metrics <- function(X,
                           bcbc_run,
                           weighted_clusters = TRUE,
                           percent_noise = 0.25,
                           num_row_clusters = NA,
                           num_col_clusters = NA,
                           rss_offset = sd(X) / 100) {
  U <- bcbc_run$U
  w <- bcbc_run$w
  lambda <- bcbc_run$lambda
  n <- nrow(X)
  p <- ncol(X)
  swept <- sweep(U, 2, w ^ 2 + lambda * w, "*")
  rss <- sum(sweep(X - U, 2, w ^ 2 + lambda * w, "*") ^ 2) + rss_offset
  rss_no_lambda_sq <- sum(sweep(X - U, 2, w ^ 2, "*") ^ 2) + rss_offset
  rss_no_lambda_no_sq <- sum(sweep(X - U, 2, w, "*") ^ 2) + rss_offset
  if (weighted_clusters) {
    solution_mat <- swept
  } else {
    solution_mat <- U
  }

  if(is.na(num_row_clusters)) {
    row_threshold <- sd(dist(solution_mat)) * percent_noise
  } else {
    row_threshold <- get_threshold_for_k_components(dist(solution_mat),
                                                    num_row_clusters)
  }
  if(is.na(num_col_clusters)) {
    col_threshold <- sd(dist(t(solution_mat))) * percent_noise
  } else {
    col_threshold <- get_threshold_for_k_components(dist(t(solution_mat)),
                                                    num_col_clusters)
  }
  row_clusters <- get_row_clusters(solution_mat, row_threshold)
  col_clusters <- get_row_clusters(t(solution_mat), col_threshold)

  num_row_clusters <- length(unique(row_clusters))
  num_col_clusters <- length(unique(col_clusters))

  sparsity_kurtosis <- sum(w ^ 2) ^ 2 / sum(w ^ 4)
  sparsity_non_zero_w <- sum(w > 1e-10)
  sparsity_non_zero_sol <- sum(abs(swept) > 1e-10)
  return(list(
    cv_metrics = list(
      num_row_clusters = num_row_clusters,
      num_col_clusters = num_col_clusters,
      sparsity_kurtosis = sparsity_kurtosis,
      sparsity_non_zero_w = sparsity_non_zero_w,
      sparsity_non_zero_sol = sparsity_non_zero_sol,
      rss = rss,
      rss_no_lambda_sq = rss_no_lambda_sq,
      rss_no_lambda_no_sq = rss_no_lambda_no_sq,
      BIC = n * p * log(rss_no_lambda_sq / n * p) +
        1 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC1 = n * p * log(rss_no_lambda_sq / n * p) +
        2 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC2 = n * p * log(rss_no_lambda_sq / n * p) +
        3 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      BIC_no_sq = n * p * log(rss_no_lambda_no_sq / n * p) +
        1 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC1_no_sq = n * p * log(rss_no_lambda_no_sq / n * p) +
        2 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC2_no_sq = n * p * log(rss_no_lambda_no_sq / n * p) +
        3 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w)
    ),
    row_clusters = row_clusters,
    col_clusters = col_clusters
  ))
}


#' Title
#'
#' @param X
#' @param bcbc_runs
#' @param all_params
#' @param weighted_clusters
#' @param percent_noise
#' @param num_row_clusters
#' @param num_col_clusters
#'
#' @return
#' @export
#'
#' @examples
get_all_cv_metrics <- function(X,
                               bcbc_runs,
                               all_params = data.frame(),
                               weighted_clusters = TRUE,
                               percent_noise = 0.25,
                               num_row_clusters = NA,
                               num_col_clusters = NA) {
  cv_data <- data.frame(all_params) |>
    mutate(index = row_number())
  row_clusters <- list()
  col_clusters <- list()

  for (param_set in 1:length(bcbc_runs)) {
    cv_info <- get_cv_metrics(
      X,
      bcbc_runs[[param_set]],
      weighted_clusters = weighted_clusters,
      percent_noise = percent_noise,
      num_row_clusters = num_row_clusters,
      num_col_clusters = num_col_clusters
    )

    cv_data[param_set, names(cv_info$cv_metrics)] <- cv_info$cv_metrics
    row_clusters[[param_set]] = cv_info$row_clusters
    col_clusters[[param_set]] = cv_info$col_clusters
  }

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

#' Title
#'
#' @param X
#' @param lambdas
#' @param k_rows
#' @param k_cols
#' @param gammas
#' @param tmaxs
#' @param tmax_outers
#' @param tmax_cobras
#' @param phis
#' @param recalculate_weights
#' @param tols
#' @param weighted_clusters
#' @param percent_noise
#' @param progress
#' @param ... passed to BCBC
#'
#' @return
#' @import future.apply
#' @export
#'
#' @examples
cv.BCBC <- function(X,
                    lambdas = c(1),
                    k_rows = c(2),
                    k_cols = c(2),
                    gammas = c(1),
                    tmaxs = c(NA),
                    tmax_outers = c(100),
                    tmax_cobras = c(100),
                    phis = c(0.5),
                    recalculate_weights = c(TRUE),
                    tols = c(1e-6),
                    weighted_clusters = TRUE,
                    percent_noise = c(0.1),
                    num_row_clusters = NA,
                    num_col_clusters = NA,
                    ...) {

  all_params <-
    expand.grid(
      lambda = lambdas,
      k_row = k_rows,
      k_col = k_cols,
      gamma = gammas,
      tmax = tmaxs,
      tmax_outer = tmax_outers,
      tmax_cobra = tmax_cobras,
      phi = phis,
      recalculate_weights = recalculate_weights,
      tol = tols
    )

  cv_data <- data.frame(all_params) |>
    mutate(index=row_number())

  all_runs <- list()
  row_clusters <- list()
  col_clusters <- list()

  future_results <- future_lapply(
    1:nrow(all_params),
    future.seed = TRUE,
    function(param_set) {
      params <- all_params[param_set,]
      result <- do.call(BCBC,
                        c(list(X=X), params, ...))

      cv_metrics <- get_cv_metrics(X,
                                   result,
                                   weighted_clusters = weighted_clusters,
                                   percent_noise = percent_noise,
                                   num_row_clusters = num_row_clusters,
                                   num_col_clusters = num_col_clusters)

      list(
        result = result,
        cv_metrics = cv_metrics$cv_metrics,
        row_clusters = cv_metrics$row_clusters,
        col_clusters = cv_metrics$col_clusters
      )
    })

  for(param_set in 1:nrow(all_params)) {
    all_runs[[param_set]] = future_results[[param_set]]$result
    cv_metrics <- future_results[[param_set]]$cv_metrics
    cv_data[param_set, names(cv_metrics)] <- cv_metrics
    row_clusters[[param_set]] = future_results[[param_set]]$row_clusters
    col_clusters[[param_set]] = future_results[[param_set]]$col_clusters
  }

  return(list(all_runs = all_runs, cv_data = cv_data,
              row_clusters = row_clusters, col_clusters = col_clusters,
              all_params = all_params))
}

melt_cv_data <- function(tuning_data) {
  # Makes a big df mostly for visualization
  num_runs <- nrow(tuning_data$cv_data)
  to_return <-
    cbind(matrix_fit_to_df(tuning_data$all_runs[[1]]), tuning_data$cv_data[1, ])
  for (i in 2:num_runs) {
    to_return <- rbind(to_return,
                       cbind(matrix_fit_to_df(tuning_data$all_runs[[i]]), tuning_data$cv_data[i, ]))
  }
  return(to_return)
}
