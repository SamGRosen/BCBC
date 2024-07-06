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
#' @param return_edges
#'
#' @return
#' @importFrom dbscan kNN
#' @import Matrix
#' @export
#'
#' @examples
fast_gkn_weights <- function(X,
                             k_row,
                             k_col,
                             phi,
                             approx = 0,
                             return_edges = TRUE) {
  p <- ncol(X)
  n <- nrow(X)

  all_row_knn <- kNN(X, k_row, sort = FALSE, approx = approx)
  all_col_knn <- kNN(t(X), k_col, sort = FALSE, approx = approx)

  # all_row_knn <- FNN::get.knn(X, k_row, algorithm="cover_tree")
  # all_col_knn <- FNN::get.knn(t(X), k_col, algorithm="cover_tree")
  unique_row_edges <- data.frame(all_row_knn$id) |>
    mutate(index = row_number()) |>
    pivot_longer(!index) |>
    mutate(
      neighbor_index = as.integer(str_remove(name, "X")),
      from_node = pmin(index, value),
      to_node = pmax(index, value)
    ) |> # TODO: May be worth double counting some edges
    distinct(from_node, to_node, .keep_all = TRUE)

  unique_col_edges <- data.frame(all_col_knn$id) |>
    mutate(index = row_number()) |>
    pivot_longer(!index) |>
    mutate(
      neighbor_index = as.integer(str_remove(name, "X")),
      from_node = pmin(index, value),
      to_node = pmax(index, value)
    ) |>
    distinct(from_node, to_node, .keep_all = TRUE)

  all_row_dist_sq <- all_row_knn$dist[cbind(unique_row_edges$from_node,
                                            unique_row_edges$neighbor_index)]^2
  all_col_dist_sq <- all_col_knn$dist[cbind(unique_col_edges$from_node,
                                            unique_col_edges$neighbor_index)]^2

  row_weights <- exp(-all_row_dist_sq * phi / p)
  row_weights <- row_weights / sum(row_weights) / sqrt(p)
  row_fusion <- sum(row_weights * sqrt(all_row_dist_sq))
  col_weights <- exp(-all_col_dist_sq * phi / n)
  col_weights <- col_weights / sum(col_weights) / sqrt(n)
  col_fusion <- sum(col_weights * sqrt(all_col_dist_sq))

  # Construct edge-incidence matrices
  E_row <-
    create_edge_incidence_edges(cbind(unique_row_edges$from_node, unique_row_edges$to_node),
                                n)
  E_col <-
    create_edge_incidence_edges(cbind(unique_col_edges$from_node, unique_col_edges$to_node),
                                p)
  to_return <- list(
    w_row = row_weights,
    w_col = col_weights,
    E_row = E_row,
    E_col = E_col,
    row_fusion = row_fusion,
    col_fusion = col_fusion
  )

  if(return_edges) {
    to_return$unique_row_edges = unique_row_edges
    to_return$unique_col_edges = unique_col_edges
  }

  return(to_return)
}
