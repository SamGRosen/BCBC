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


#' Title
#'
#' @param bcbc_result
#' @param percent_of_noise
#' @param cluster_w_weights
#'
#' @return
#' @export
#'
#' @examples
thresholded_solution <- function(bcbc_result,
                                 percent_of_noise,
                                 cluster_w_weights = TRUE) {
  to_return <- list()
  w <- bcbc_result$w
  U <- bcbc_result$U
  lambda <- bcbc_result$lambda
  if (cluster_w_weights) {
    w_vals <- sqrt(w + lambda * w)
    solution_mat <- sweep(U, 2, w_vals, "*")
  } else {
    solution_mat <- U
  }
  row_sd <- sd(sqrt(rowSums(solution_mat ^ 2)))
  threshold <- row_sd * percent_of_noise
  clustering <- centroid_rows(U, solution_mat, threshold)
  bcbc_result$U <- clustering$mat
  bcbc_result$cluster_info <- clustering$cluster_info

  to_return$w <- w
  to_return$lambda <- lambda
  to_return$cluster_info <- clustering$cluster_info
  to_return$U <- clustering$mat

  return(to_return)
}
