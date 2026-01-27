#' Build row and column info for usage with `BCBC`, `cobra` and subroutines
#' using specification from paper. Use `knn_graph` or `knn_graph_approx` for
#' more precise control.
#'
#' @param X data matrix
#' @param k_row number of nearest neighbors for rows
#' @param k_col number of nearest neighbors for columns
#' @param approximate use RcppHNSW for calculation
#' @param hnsw_args passed to `RcppHNSW::hnsw_knn` if `approximate` is true
#' @seealso [knn_graph()]
#' @seealso [knn_graph_approx()]
#' @return list of row and column info similar to `knn_graph` return info under
#' `W_row` and `W_col` keys.
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(n = 100,
#'                             p = 100,
#'                             num_row_clusters = 5,
#'                             num_col_clusters = 5,
#'                             p_extra = 25,
#'                             shuffle = FALSE)
#' fusion_wts <- fast_gkn_weights(t(checker$X), k_row = 4, k_col = 4)
#' calc_fusion_term(fusion_wts$unique_row_edges, fusion_wts$w_row, checker$X, rows=FALSE)
#' calc_fusion_term(fusion_wts$unique_col_edges, fusion_wts$w_col, checker$X, rows=TRUE)
fast_gkn_weights <- function(X,
                             k_row,
                             k_col,
                             approximate = FALSE,
                             hnsw_args = list(ef = 50)) {
  if (approximate) {
    row_knn <- knn_graph_approx(
      X,
      k_row,
      sum_to_one = TRUE,
      rescale_by_dim = ncol(X),
      uniform = FALSE,
      phi = 1 / ncol(X),
      hnsw_args = hnsw_args
    )
    col_knn <- knn_graph_approx(
      t(X),
      k_col,
      sum_to_one = TRUE,
      rescale_by_dim = nrow(X),
      uniform = FALSE,
      phi = 1 / nrow(X),
      hnsw_args = hnsw_args
    )
  } else {
    row_knn <- knn_graph(
      X,
      k_row,
      sum_to_one = TRUE,
      rescale_by_dim = ncol(X),
      uniform = FALSE,
      phi = 1 / ncol(X)
    )
    col_knn <- knn_graph(
      t(X),
      k_col,
      sum_to_one = TRUE,
      rescale_by_dim = nrow(X),
      uniform = FALSE,
      phi = 1 / nrow(X)
    )
  }

  list(
    W_row = row_knn,
    W_col = col_knn
  )
}

#' Calculate exact nearest neighbors for use with `BCBC`, `cobra` and subroutines.
#'
#' @param X data matrix
#' @param k number of nearnest neighbors
#' @inheritParams build_weights_for_edges
#' @inherit build_weights_for_edges return
#'
#' @importFrom dbscan kNN
#' @export
#'
knn_graph <- function(X,
                      k,
                      sum_to_one = TRUE,
                      rescale_by_dim = ncol(X),
                      uniform = FALSE,
                      phi = 1 / ncol(X)) {
  all_row_knn <- dbscan::kNN(X, k, sort = FALSE)

  build_weights_for_edges(
    all_row_knn$id,
    all_row_knn$dist,
    sum_to_one = sum_to_one,
    rescale_by_dim = rescale_by_dim,
    uniform = uniform,
    phi = phi
  )
}


#' Calculate approximate nearest neighbors for use with `BCBC`, `cobra` and subroutines.
#'
#' @param X data matrix
#' @param k number of nearest neighbors
#' @param hnsw_args passed to `RcppHNSW::hnsw_knn`
#' @inheritParams build_weights_for_edges
#' @inherit build_weights_for_edges return
#' @import RcppHNSW
#' @export
#'
knn_graph_approx <- function(X,
                             k,
                             sum_to_one = TRUE,
                             rescale_by_dim = ncol(X),
                             uniform = FALSE,
                             phi = 1 / ncol(X),
                             hnsw_args = list()) {
  knn_results <- do.call(hnsw_knn, c(list(
    X = X, k = k + 1, byrow = TRUE  # k + 1 to avoid counting self-loops
  ), hnsw_args))

  build_weights_for_edges(
    knn_results$idx,
    knn_results$dist,
    sum_to_one = sum_to_one,
    rescale_by_dim = rescale_by_dim,
    uniform = uniform,
    phi = phi
  )
}


#' Build completely connected affinity graph
#'
#' @param X data matrix
#' @param sum_to_one sum the weights to one
#' @param rescale_by_dim divide the weights by sqrt(ncol(X))
#' @param uniform keep all weights uniform
#' @param phi exponential kernel scaling
#'
#' @return sparseweights object
#' @export
#'
full_connectivity <- function(X,
                              sum_to_one = TRUE,
                              rescale_by_dim = TRUE,
                              uniform = FALSE,
                              phi = 1 / ncol(X)) {
  n <- nrow(X)
  p <- ncol(X)

  # Left column should be greater than right for CCMMR
  all_edges <- t(combn(n:1, 2))
  all_dist <- c(dist(X))

  weights <- distances_to_weights(all_dist,
                                  sum_to_one = sum_to_one,
                                  rescale_by_dim = p,
                                  uniform = uniform,
                                  phi = phi)

  structure(
    list(values = weights, keys = all_edges),
    class = "sparseweights"
  )
}

#' Build list for usage with `cobra`.
#'
#' @param edges id attribute matrix output from `dbscan::kNN` or (idx) `RcppHNSW::hnsw_knn`
#' @param distances dist attribute output from `dbscan::kNN` or `RcppHNSW::hnsw_knn`
#' @param uniform treat all edges as equal magnitude
#' @param rescale_by_dim rescale weights by dividing by sqrt(this number) (if finite)
#' @param phi bandwidth for Gaussian kernel
#' @param sum_to_one make weights sum to one before any dimension rescaling
#'
#' @return `sparseweights` object (list) with
#'   1. `keys` matrix of edges
#'   1. `values` weights for edges in `keys`
#'   1. `fusion` sum of fusion term when using the weights and distances
#' @export
build_weights_for_edges <- function(edges,
                                    distances,
                                    sum_to_one = TRUE,
                                    rescale_by_dim = NA,
                                    uniform = FALSE,
                                    phi = 1) {
  unique_edges <- data.frame(edges) |>
    mutate(point_index = row_number()) |>
    pivot_longer(!point_index, names_to = "kth_neighbor", values_to = "neighbor_index") |>
    mutate(
      kth_neighbor = as.integer(substring(kth_neighbor, 2)),
      from_node = pmax(point_index, neighbor_index),
      to_node = pmin(point_index, neighbor_index)
    ) |>
    filter(from_node != to_node) |>  # Remove self loops
    distinct(from_node, to_node, .keep_all = TRUE)

  all_row_dist <- distances[cbind(unique_edges$from_node,
                                     unique_edges$kth_neighbor)]

  keys = matrix(c(unique_edges$from_node, unique_edges$to_node),
                nrow = nrow(unique_edges), ncol = 2)

  weights = distances_to_weights(
    all_row_dist,
    sum_to_one = sum_to_one,
    rescale_by_dim = rescale_by_dim,
    uniform = uniform,
    phi = phi
  )

  fusion <- sum(weights * all_row_dist)

  structure(
    list(values = weights, keys = keys, fusion = fusion),
    class = "sparseweights"
  )
}

distances_to_weights <- function(distances,
                                 sum_to_one = TRUE,
                                 rescale_by_dim = NA,
                                 uniform = FALSE,
                                 phi = 1) {
  if(!uniform) {
    weights <- exp(-phi * distances^2)
  } else {
    weights <- rep(1, length(distances))
  }

  if(sum_to_one) {
    weights <- weights / sum(weights)
  }

  if(is.finite(rescale_by_dim) && rescale_by_dim >= 1) {
    weights <- weights / sqrt(rescale_by_dim)
  }

  weights
}
