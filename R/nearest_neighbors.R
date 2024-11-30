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

#' Title
#'
#' @param X
#' @param k_row
#' @param k_col
#' @param phi
#' @param approx
#' @param rescale
#' @param uniform
#' @param hnsw_args
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param k
#' @param rescale
#' @param uniform
#' @param phi
#'
#' @return
#' @importFrom dbscan kNN
#' @export
#'
#' @examples
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


#' Title
#'
#' @param X
#' @param k
#' @param rescale
#' @param uniform
#' @param phi
#' @param byrow
#' @param hnsw_args
#'
#' @return
#' @import RcppHNSW
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param rescale
#' @param uniform
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param X
#' @param rescale
#' @param uniform
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param X
#' @param cluster_assignments
#'
#' @return
#' @export
#'
#' @examples
connected_oracle <- function(X, cluster_assignments) {
  n <- nrow(X)
  cluster_labels <- unique(cluster_assignments)

  unique_edges <- tibble(
    from_node = numeric(),
    to_node = numeric(),
    distance = numeric()
  )
  for(cluster in cluster_labels) {
    cluster_indices <- which(cluster_assignments == cluster)
    inside_cluster_edges <- combn(cluster_indices, 2)
    unique_edges <- unique_edges |>
      rbind(tibble(
        from_node = inside_cluster_edges[1, ],
        to_node = inside_cluster_edges[2, ],
        distance = 0
      ))
  }
  first_cluster_indices <- rep(0, length(cluster_labels))
  i <- 1
  for(cluster in cluster_labels) {
    first_cluster_indices[i] <- min(which(cluster_assignments == cluster))
    i <- i + 1
  }

  bridges <- combn(first_cluster_indices, 2)
  for(bridge_index in 1:ncol(bridges)) {
    unique_edges <- unique_edges |>
      add_row(from_node = bridges[1, bridge_index],
              to_node = bridges[2, bridge_index],
              distance = sqrt(sum((X[bridges[1, bridge_index], ] -
                                     X[bridges[2, bridge_index], ])^2)))
  }

  # Construct edge-incidence matrices
  E_row <- create_edge_incidence_edges(
    cbind(unique_edges$from_node, unique_edges$to_node),
    n)

  list(
    weights = rep(1, nrow(unique_edges)),
    E_row = E_row,
    fusion = sum(unique_edges$distance),
    unique_edges = unique_edges
  )
}


#' Title
#'
#' @param edges
#' @param distances
#' @param uniform
#' @param rescale
#' @param phi
#'
#' @return
#' @export
#'
#' @examples
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
