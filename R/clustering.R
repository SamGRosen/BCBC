#' Get row clusters memberships with thresholding radius.
#'
#' @param mat to find row memberships for
#' @param threshold radius to join rows into communities
#'
#' @return vector of row memberships
#' @importFrom dbscan comps
#' @importFrom dbscan frNN
#' @export
#' @seealso [dbscan::frNN()]
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 100, 5, shuffle = TRUE)
#' get_row_clusters(checker$X, 10*sd(checker$X)^2)
get_row_clusters <- function(mat, threshold) {
  if (is.null(threshold) ||
      is.na(threshold) || threshold < 0) {
    warning(paste("threshold is invalid:", threshold))
    threshold <- 0
  }
  row_clusters <- dbscan::comps(dbscan::frNN(mat, threshold, sort=F))
  names(row_clusters) <- 1:nrow(mat)

  row_clusters
}


#' Get radius to create k connected components in a nearest neighbor graph
#'
#' @param dist_mat distance matrix of data points
#' @param k number of connected components
#' @param max_threshold highest radius to use
#' @param maxiter number of iters used to find radius
#'
#' @return radius for k connected components
#' @importFrom dbscan comps
#' @importFrom dbscan frNN
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 0, shuffle = TRUE)
#' threshold_for_k_components(dist(checker$X), 20)
threshold_for_k_components <- function(dist_mat,
                                       k,
                                       max_threshold = max(dist_mat),
                                       maxiter = 50) {
  get_num_components_err <- function(threshold) {
    length(unique(dbscan::comps(dbscan::frNN(
      dist_mat, threshold, sort = F
    )))) - k
  }

  if (get_num_components_err(max_threshold) >= 0) {
    return(max_threshold)
  }

  root <- uniroot(get_num_components_err,
                  c(0, max_threshold),
                  maxiter = maxiter,
                  tol = 1e-10)

  root$root
}

#' Calculate biclusters according to some linking criteria
#'
#' @param U matrix to bicluster
#' @param weights of feature importance
#' @param lambda used to complete fit
#' @param percent_noise percent of distance standard deviation used as joining radius
#' @param num_row_clusters try to find this many row clusters if defined
#' @param num_col_clusters try to find this many column clusters if defined
#'
#' @return list with items
#'   1. `row_clusters` row membership index
#'   1. `col_clusters` column membership index
#'   1. `row_threshold` radius used to calculate row memberships
#'   1. `col_threshold` radius used to calculate column memberships
#' @seealso [unweighted_bicluster_assignments()]
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.05, gamma = 200)
#' assignments <- bicluster_assignments(
#'   bcbc_fit$U,
#'   weights = bcbc_fit$w,
#'   lambda = bcbc_fit$lambda,
#'   percent_noise = 0.1
#' )
#' as_df <- bcbc_result_to_df(bcbc_fit, labels = assignments$row_clusters)
#' plot_fit(as_df)
bicluster_assignments <- function(U,
                                  weights=NA,
                                  lambda=NA,
                                  percent_noise=0.25,
                                  num_row_clusters=NA,
                                  num_col_clusters=NA) {
  n = nrow(U)
  p = ncol(U)

  non_zero <- which(weights > 0)
  if(length(non_zero) < 2) {
    warning(paste("Given weights are invalid", paste(weights, collapse = " ")))
    weights = rep(0, p)
    weights[1] = 1/2
    weights[2] = 1/2
    non_zero <- which(weights > 0)
  }

  if(any(is.na(weights))) {
    weights <- rep(1, p)
    lambda <- NA
  }

  non_zero_w <- weights[non_zero]

  if(is.finite(lambda)) {
    non_zero_w <- non_zero_w^2 + lambda * non_zero_w
  }

  removed <- U[, non_zero]
  weighted_and_removed <- sweep(removed, 2, sqrt(non_zero_w), "*")

  if(is.na(num_row_clusters)) {
    row_threshold <- sd(dist(weighted_and_removed)) * percent_noise
  } else {
    row_threshold <- threshold_for_k_components(dist(weighted_and_removed),
                                                    num_row_clusters)
  }
  if(is.na(num_col_clusters)) {
    col_threshold <- sd(dist(t(removed))) * percent_noise
  } else {
    col_threshold <- threshold_for_k_components(dist(t(removed)),
                                                    num_col_clusters)
  }
  row_clusters <- get_row_clusters(weighted_and_removed, row_threshold)
  col_clusters_non_zero <- get_row_clusters(t(removed), col_threshold)

  col_clusters <- rep(0, p)  # 0 is dummy label
  col_clusters[non_zero] <- col_clusters_non_zero

  list(
    row_clusters = row_clusters,
    col_clusters = col_clusters,
    row_threshold = row_threshold,
    col_threshold = col_threshold
  )
}

#' Calculate biclusters according to some linking criteria without feature weights
#'
#' @inheritParams bicluster_assignments
#' @inherit bicluster_assignments return
#' @seealso [bicluster_assignments()]
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.05, gamma = 200)
#' assignments <- unweighted_bicluster_assignments(
#'   bcbc_fit$U,
#'   percent_noise = 0.1
#' )
#' as_df <- bcbc_result_to_df(bcbc_fit, labels = assignments$row_clusters)
#' plot_fit(as_df)
unweighted_bicluster_assignments <- function(U,
                                             num_row_clusters = NA,
                                             num_col_clusters = NA,
                                             percent_noise = 0.25) {

  n = nrow(U)
  p = ncol(U)

  if(is.na(num_row_clusters)) {
    row_threshold <- sd(dist(U)) * percent_noise
  } else {
    row_threshold <- threshold_for_k_components(dist(U),
                                                    num_row_clusters)
  }
  if(is.na(num_col_clusters)) {
    col_threshold <- sd(dist(t(U))) * percent_noise
  } else {
    col_threshold <- threshold_for_k_components(dist(t(U)),
                                                    num_col_clusters)
  }
  row_clusters <- get_row_clusters(U, row_threshold)
  col_clusters_non_zero <- get_row_clusters(t(U), col_threshold)

  col_clusters <- col_clusters_non_zero

  list(
    row_clusters = row_clusters,
    col_clusters = col_clusters,
    row_threshold = row_threshold,
    col_threshold = col_threshold
  )
}

#' Calculate least square estimate for bicluster centers given assignments
#'
#' @param X data matrix to find centers for
#' @param row_clusters cluster membership of rows
#' @param col_clusters cluster membership of columns
#'
#' @return list with entries `fit` for fitted centers and `assignments` for bicluster labels
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = FALSE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.05, gamma = 200)
#' assignments <- bicluster_assignments(
#'   bcbc_fit$U,
#'   weights = bcbc_fit$w,
#'   lambda = bcbc_fit$lambda,
#'   percent_noise = 0.1
#' )
#' centers <- bicluster_centers(
#'   checker$X,
#'   assignments$row_clusters,
#'   assignments$col_clusters
#' )
#' plot_fit(matrix_fit_to_df(centers$fit))
bicluster_centers <- function(X, row_clusters, col_clusters) {
  bicluster_assignments <- outer(row_clusters, col_clusters, paste)
  bicluster_assignments[, col_clusters == 0] <- "0, 0"

  as_df <- tibble(
    row = rep(1:nrow(X), each = ncol(X)),
    col = rep(1:ncol(X), nrow(X)),
    assignment = as.vector(bicluster_assignments),
    value = as.vector(X)
  )

  num_groups <- length(unique(as_df$assignment))
  group_means <- as_df |>
    group_by(assignment) |>
    mutate(val = mean(value)) |>
    ungroup()

  fitted <- matrix(group_means$val, nrow=nrow(X), ncol=ncol(X))
  return(list(
    fit = fitted,
    assignments = bicluster_assignments
  ))
}
