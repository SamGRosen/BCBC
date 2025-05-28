create_edge_incidence_edges <- function(P, n) {
  nEdges <- nrow(P)
  E <- Matrix(
    data = 0,
    nrow = nEdges,
    ncol = n,
    sparse = TRUE
  )
  r <- 1:nEdges
  col <- P[, 1]
  E[(col - 1) * nEdges + r] <- 1
  col <- P[, 2]
  E[(col - 1) * nEdges + r] <- -1
  return(E)
}

#' Build row and column info for usage with `BCBC`, `cobra` and subroutines.
#'
#' @param X data matrix
#' @param k_row number of nearest neighbors for rows
#' @param k_col number of nearest neighbors for columns
#' @param approximate use RcppHNSW for calculation
#' @param hnsw_args passed to `RcppHNSW::hnsw_knn`
#' @inheritParams build_weights_for_edges
#' @seealso [knn_graph()]
#' @seealso [knn_graph_approx()]
#' @return list of row and column info similar to `knn_graph` return info
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(n = 100,
#'                             p = 100,
#'                             num_row_clusters = 5,
#'                             num_col_clusters = 5,
#'                             p_extra = 25,
#'                             shuffle = FALSE)
#' fusion_wts <- fast_gkn_weights(t(checker$X), k_row = 4, k_col = 4, phi = 1)
#' calc_fusion_term(fusion_wts$unique_row_edges, fusion_wts$w_row, checker$X, rows=FALSE)
#' calc_fusion_term(fusion_wts$unique_col_edges, fusion_wts$w_col, checker$X, rows=TRUE)
fast_gkn_weights <- function(X, k_row, k_col, phi,
                             approximate = FALSE,
                             rescale = TRUE,
                             uniform = FALSE,
                             hnsw_args = list(ef = 50)) {
  if(approximate) {
    row_knn <- knn_graph_approx(X, k_row, rescale, uniform, phi, hnsw_args, byrow = TRUE)
    col_knn <- knn_graph_approx(X, k_row, rescale, uniform, phi, hnsw_args, byrow = FALSE)
  } else {
    row_knn <- knn_graph(X, k_row, rescale, uniform, phi)
    col_knn <- knn_graph(t(X), k_col, rescale, uniform, phi)
  }
  list(
    w_row = row_knn$weights,
    w_col = col_knn$weights,
    E_row = row_knn$E_row,
    E_col = col_knn$E_row,
    row_fusion = row_knn$fusion,
    col_fusion = col_knn$fusion,
    unique_row_edges = row_knn$unique_edges,
    unique_col_edges = col_knn$unique_edges
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
                      rescale = FALSE,
                      uniform = FALSE,
                      phi = 1) {
  all_row_knn <- kNN(X, k, sort = FALSE)
  if (is.finite(rescale) && rescale > 0) {
    rescale <- ncol(X)
  }
  build_weights_for_edges(
    all_row_knn$id,
    all_row_knn$dist,
    uniform = uniform,
    rescale = rescale,
    phi = phi
  )
}


#' Calculate approximate nearest neighbors for use with `BCBC`, `cobra` and subroutines.
#'
#' @param X data matrix
#' @param k number of nearest neighbors
#' @param byrow calculate for rows instead of columns
#' @param hnsw_args passed to `RcppHNSW::hnsw_knn`
#' @inheritParams build_weights_for_edges
#' @inherit build_weights_for_edges return
#' @import RcppHNSW
#' @export
#'
knn_graph_approx <- function(X,
                             k,
                             rescale = FALSE,
                             uniform = FALSE,
                             phi = 1,
                             byrow = FALSE,
                             hnsw_args = list()) {

  if(byrow) {
    knn_results <- do.call(hnsw_knn, c(list(
      X = t(X), k = k + 1, byrow = FALSE  # k + 1 to avoid counting self-loops
    ), hnsw_args))
  } else {
    knn_results <- do.call(hnsw_knn, c(list(
      X = X, k = k + 1, byrow = FALSE
    ), hnsw_args))
  }

  if (is.finite(rescale) && rescale > 0) {
    if(byrow) {
      rescale <- ncol(X)
    } else {
      rescale <- nrow(X)
    }
  }
  build_weights_for_edges(
    t(knn_results$idx),
    t(knn_results$dist),
    uniform = uniform,
    rescale = rescale,
    phi = phi
  )
}

#' Find geometric mediod and build star affinity graph
#'
#' @param X data matrix
#' @inheritParams build_weights_for_edges
#' @inherit build_weights_for_edges return
#'
#' @export
#'
geometric_medoid <- function(X,
                             rescale = FALSE,
                             uniform = FALSE,
                             phi = 1) {
  n <- nrow(X)
  p <- ncol(X)

  all_dist <- as.matrix(dist(X))
  medoid_index <- which.min(rowSums(all_dist))

  # Add duplicate edges as they are removed anyways
  all_edges <- matrix(medoid_index, nrow = n, ncol = n - 1)
  all_edges[medoid_index, ] <- (1:n)[-medoid_index]

  # Set to have duplicate distances, except at medoid
  all_dist_medoid <- matrix(all_dist[medoid_index,], nrow=n, ncol=n - 1)
  all_dist_medoid[medoid_index, ] <- all_dist[medoid_index, -medoid_index]

  if (is.finite(rescale) && rescale > 0) {
    rescale <- p
  }
  build_weights_for_edges(
    all_edges,
    all_dist_medoid,
    uniform = uniform,
    rescale = rescale,
    phi = phi
  )
}


#' Build completely connected affinity graph
#'
#' @param X data matrix
#' @inheritParams build_weights_for_edges
#' @inherit build_weights_for_edges return
#' @export
#'
full_connectivity <- function(X,
                              rescale = NA,
                              uniform = FALSE,
                              phi = 1) {
  n <- nrow(X)
  p <- ncol(X)

  all_edges <- matrix(1:nrow(X),
                      nrow = n,
                      ncol = n,
                      byrow = TRUE)
  all_dist <- as.matrix(dist(X))

  # Remove loops https://stackoverflow.com/a/18879755
  diag(all_edges) <- NA
  all_edges <- t(matrix(t(all_edges)[which(!is.na(all_edges))],
                        nrow = n - 1, ncol = n))

  diag(all_dist) <- NA
  all_dist <- t(matrix(t(all_dist)[which(!is.na(all_dist))],
                       nrow = n - 1, ncol = n))


  if (is.finite(rescale) && rescale > 0) {
    rescale <- p
  }
  build_weights_for_edges(
    all_edges,
    all_dist,
    uniform = uniform,
    rescale = rescale,
    phi = phi
  )
}

#' Build list for usage with `cobra`.
#'
#' @param edges data frame describing
#' @param distances corresponding to each edge
#' @param uniform treat all edges as equal magnitude
#' @param rescale weights and bandwidth by this factor
#' @param phi bandwidth for Gaussian kernel
#'
#' @return list with
#'   1. `weights` vector of edge weights
#'   1. `E_row` sparse unweighted edge incidence matrix
#'   1. `fusion` sum of fusion term when using the weights and distances
#'   1. `unique_edges` data frame describing edge information
#' @export
build_weights_for_edges <- function(edges,
                                    distances,
                                    uniform = FALSE,
                                    rescale = 1,
                                    phi = 1) {
  unique_edges <- data.frame(edges) |>
    mutate(index = row_number()) |>
    pivot_longer(!index) |>
    mutate(
      neighbor_index = as.integer(sub("X", "", name)),
      from_node = pmin(index, value),
      to_node = pmax(index, value)
    ) |>
    filter(from_node != to_node) |>  # Remove self loops
    distinct(from_node, to_node, .keep_all = TRUE)

  all_row_dist_sq <- distances[cbind(unique_edges$from_node,
                                     unique_edges$neighbor_index)]^2
  if (!uniform) {
    if (is.finite(rescale) && rescale > 1) {
      phi <- phi / rescale
    }
    weights <- exp(-all_row_dist_sq * phi)
  } else {
    weights <- rep(1, nrow(unique_edges))
  }

  if(is.finite(rescale) && rescale > 0) {
    weights <- weights / sum(weights) / sqrt(rescale)
  }

  fusion <- sum(weights * sqrt(distances[cbind(unique_edges$from_node,
                                               unique_edges$neighbor_index)]^2))

  # Construct edge-incidence matrices
  E_row <-
    create_edge_incidence_edges(
      cbind(unique_edges$from_node, unique_edges$to_node),
      max(unique_edges$from_node, unique_edges$to_node)
    )

  list(
    weights = weights,
    E_row = E_row,
    fusion = fusion,
    unique_edges = unique_edges
  )
}
