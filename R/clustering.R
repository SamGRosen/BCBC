#' Title
#'
#' @param mat
#' @param threshold
#'
#' @return
#' @import dbscan
#' @export
#'
#' @examples
get_row_clusters <- function(mat, threshold) {
  if (is.null(threshold) ||
      is.na(threshold) || threshold < 0) {
    warning(paste("threshold is invalid:", threshold))
    threshold <- 0
  }
  row_clusters <- dbscan::comps(frNN(mat, threshold, sort=F))
  names(row_clusters) <- 1:nrow(mat)

  row_clusters
}


#' Get radius to create k connected components in a nearest neighbor graph
#'
#' @param dist_mat
#' @param k
#' @param max_threshold
#' @param maxiter
#'
#' @return
#' @import dbscan
#' @export
#'
#' @examples
threshold_for_k_components <- function(dist_mat,
                                       k,
                                       max_threshold = max(dist_mat),
                                       maxiter = 50) {
  get_num_components_err <- function(threshold) {
    length(unique(dbscan::comps(frNN(
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

#' Title
#'
#' @param U
#' @param weights
#' @param lambda
#' @param num_row_clusters
#' @param num_col_clusters
#' @param percent_noise
#'
#' @return
#' @export
#'
#' @examples
bicluster_assignments <- function(U,
                                  weights=NA,
                                  lambda=NA,
                                  num_row_clusters=NA,
                                  num_col_clusters=NA,
                                  percent_noise=0.25) {
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

#' Title
#'
#' @param U
#' @param num_row_clusters
#' @param num_col_clusters
#' @param percent_noise
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param row_clusters
#' @param col_clusters
#'
#' @return
#' @export
#'
#' @examples
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
