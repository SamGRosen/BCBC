centroid_rows <- function(mat, mat_for_dist, threshold, calculate_centroids=TRUE) {
  to_return <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  row_clusters <- dbscan::comps(frNN(mat_for_dist, threshold, sort=F))
  names(row_clusters) <- 1:nrow(mat)
  sorted_membership <- sort(row_clusters)
  node_indices <- as.integer(names(sorted_membership))

  if(!calculate_centroids) {
    return(list(mat = to_return, cluster_info = row_clusters))
  }

  cluster_sizes <- table(row_clusters)
  curr_index <- 1

  for (cluster_id in names(cluster_sizes)) {
    cluster_size <- cluster_sizes[cluster_id]
    cluster_members <-
      node_indices[curr_index:(curr_index + cluster_size - 1)]
    if (cluster_size > 1) {
      centroid <- colSums(mat[cluster_members,]) / cluster_size
      to_return[cluster_members,] <-
        rep(centroid, each = cluster_size)
    }

    curr_index <- curr_index + cluster_size
  }

  return(list(mat = to_return, cluster_info = row_clusters))
}


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
  clustering <- centroid_rows(mat, mat, threshold, calculate_centroids = FALSE)

  return(clustering$cluster_info)
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
get_threshold_for_k_components <- function(dist_mat,
                                           k,
                                           max_threshold = sd(dist_mat),
                                           maxiter = 25) {
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
get_weighted_biclusters <- function(U,
                                    weights=NA,
                                    lambda=NA,
                                    num_row_clusters=NA,
                                    num_col_clusters=NA,
                                    percent_noise=0.25) {
  n = nrow(U)
  p = ncol(U)
  if(any(is.na(weights))) {
    weights <- rep(1, p)
    lambda <- NA
  }
  non_zero <- which(weights > 0)
  non_zero_w <- weights[non_zero]

  if(is.finite(lambda)) {
    non_zero_w <- non_zero_w^2 + lambda * non_zero_w
  }

  removed <- U[, non_zero]
  weighted_and_removed <- sweep(removed, 2, non_zero_w, "*")

  if(is.na(num_row_clusters)) {
    row_threshold <- sd(dist(weighted_and_removed)) * percent_noise
  } else {
    row_threshold <- get_threshold_for_k_components(dist(weighted_and_removed),
                                                    num_row_clusters)
  }
  if(is.na(num_col_clusters)) {
    col_threshold <- sd(dist(t(removed))) * percent_noise
  } else {
    col_threshold <- get_threshold_for_k_components(dist(t(removed)),
                                                    num_col_clusters)
  }
  row_clusters <- get_row_clusters(weighted_and_removed, row_threshold)
  col_clusters_non_zero <- get_row_clusters(t(removed), col_threshold)

  col_clusters <- rep(0, p)  # 0 is dummy label
  col_clusters[non_zero] <- col_clusters_non_zero

  list(
    row_clusters = row_clusters,
    col_clusters = col_clusters
  )
}
