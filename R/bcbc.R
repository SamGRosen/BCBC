#' Orthogonal projection of a vector onto the simplex
#'
#' From https://github.com/kharchenkolab/vrnmf/tree/main
#' @param unproj input vector
#' @param bound sum of projected vector
#'
#' @return orthogonal projection of `unproj`
#' @export
#'
#' @examples
#' projection_onto_simplex(c(0.3, 0.3, 0.2))
#'
projection_onto_simplex <- function(unproj, bound = 1) {
  q <- sort(unproj, decreasing = TRUE, method = "quick")
  qcum <- cumsum(q)
  mu <- (qcum - bound) / 1:length(qcum)
  cond1 <- (mu[-length(mu)] - q[-1]) > 0
  if (max(cond1) == 0) {
    ind <- length(mu)
  } else {
    ind <- which.max(cond1)
  }
  return(pmax(0, unproj - mu[ind]))
}


#' Calculate the fusion term in the BCBC objective function
#'
#' @param edges `sparseweights` object
#' @param U input to fusion function
#' @param rows logical value indicating to calculate over rows
#'
#' @seealso [fast_gkn_edges()] for Gaussian k-nearest neighbor affinities
#' @seealso [knn_graph()] for Gaussian k-nearest neighbor affinities
#' @return fusion value
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(n = 100,
#'                             p = 100,
#'                             num_row_clusters = 5,
#'                             num_col_clusters = 5,
#'                             p_extra = 25,
#'                             shuffle = FALSE)
#' fusion_wts <- fast_gkn_weights(checker$X, k_row = 4, k_col = 4)
#' calc_fusion_term(fusion_wts$W_row, checker$X, rows=FALSE)
#' calc_fusion_term(fusion_wts$W_col, checker$X, rows=TRUE)
calc_fusion_term <- function(edges, U, rows=TRUE) {
  from_node <- edges$keys[, 1]
  to_node <- edges$keys[, 2]

  all_distances <- rep(NA, length(from_node))
  for(i in seq_along(from_node)) {
    if(rows) {
      all_distances[i] = sqrt(sum((U[from_node[i], ] - U[to_node[i], ])^2))
    } else {
      all_distances[i] = sqrt(sum((U[, from_node[i]] - U[, to_node[i]])^2))
    }
  }
  return(sum(edges$values * all_distances))
}


#' Perform block coordinate optimization on weights after finding fitted matrix
#'
#' @param col_sum_sq Sum of squared errors for columns
#' @param lambda hyperparameter in biconvex objective for sparsity
#'
#' @return argmin_{w in simplex} sum (w_i^2 + lambda * w_i) col_sum_sq_i
#' @export
w_coordinate_descent <- function(col_sum_sq, lambda) {
  f = function(alpha) {
    s = sum(pmax(alpha / col_sum_sq - lambda, 0))
    return(s / 2 - 1)
  }
  alpha = uniroot(f, c(0.0, 1000.0), tol = 1e-10)$root

  pmax(alpha / col_sum_sq - lambda, 0) / 2
}


#' Perform Biconvex Biclustering
#'
#' @param X data matrix
#' @param lambda sparsity hyperparameter
#' @param k_features number of nearest neighbors used to calculate feature affinity graph
#' @param k_samples number of nearest neighbors used to calculate sample affinity graph
#' @param gamma fusion hyperparameter
#' @param tmax number of max iterations for COBRA and PALM
#' @param tmax_cobra number of max iterations for COBRA
#' @param tmax_outer number of max iterations for PALM
#' @param tol for early termination
#' @param recalculate_weights if TRUE use adaptive BCBC
#' @param greedy_terminate terminate if objective function does not decrease after iteration
#' @param approx_neighbors use approximate nearest neighbors for affinity graph
#' @param hnsw_args passed to `fast_gkn_weights` for approximate NN calculation
#' @param progress show progress
#' @param fusion_wts use predetermined affinity graph
#' @param init_U use this matrix for initial PALM iterations
#'
#' @return list for BCBC fit including
#'   1. `U` fitted matrix
#'   1. `w` fitted weights
#'   1. `lambda`, `gamma`, `k_features`, `k_samples`, `recalculate_weights`  input parameters
#'   1. `cobra_diffs` difference in `U` iterates
#'   1. `w_diffs` difference in `w` iterates
#'   1. `w_path` all `w` iterates
#'   1. `objective` value of objective function for all iterates
#'   1. `rss` weighted rss value for all iterates
#'   1. `row_fusion` row fusion penalties for all iterates
#'   1. `col_fusion` col fusion penalties for all iterates
#'   1. `time` to termination
#'   1. `col_residuals` unweighted residuals per column for all iterates
#'   1. `best_obj_w` best w at lowest objective found during iteration
#'   1. `best_obj_U` best U at lowest objective found during iteration
#' @import progress
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
BCBC <- function(X,
                 lambda,
                 k_features = 4,
                 k_samples = 4,
                 gamma = 10,
                 tmax = NA, # Set this to set both below
                 tmax_cobra = 100,
                 tmax_outer = 100,
                 tol = 1e-5,
                 recalculate_weights = TRUE,
                 greedy_terminate = FALSE,
                 approx_neighbors = FALSE,
                 hnsw_args = list(),
                 progress = TRUE,
                 fusion_wts = NA,
                 init_U = X) {

  stopifnot("gamma must be > 0" = gamma > 0)
  stopifnot("lambda must be >= 0" = lambda >= 0)

  n <- dim(X)[1]
  p <- dim(X)[2]
  U <- init_U
  w <- rep(1/p, p)

  best_obj_w <- NA
  U_for_best_obj <- NA

  if (!is.na(tmax)) {
    tmax_cobra <- tmax
    tmax_outer <- tmax
  }

  cobra_diffs <- rep(NA, tmax_outer)
  w_diffs <- rep(NA, nrow = tmax_outer)
  w_path <- matrix(NA, nrow = tmax_outer, ncol = p)
  objective_vals <- rep(NA, tmax_outer)
  rss_vals <- rep(NA, tmax_outer)
  row_fusion_vals <- rep(NA, tmax_outer)
  col_fusion_vals <- rep(NA, tmax_outer)
  col_residuals <- matrix(NA, nrow = tmax_outer, ncol = p)

  if (progress) {
    pb <- txtProgressBar(min = 0,
                         max = tmax_outer,
                         initial = 0)
  }

  start <- Sys.time()

  if(!is.list(fusion_wts)) {  # If not given first iteration weights, calculate them
    fusion_wts <- fast_gkn_weights(
      X,
      k_row = k_samples,
      k_col = k_features,
      approximate = approx_neighbors,
      hnsw_args = hnsw_args
    )
  }

  for (t in 1:tmax_outer) {
    U_step <- U - sweep(U - X, 2, w ^ 2 + lambda * w, "*")

    if (recalculate_weights && t > 1) {
      fusion_wts <- fast_gkn_weights(
        U_step,  # After first step, use U_step here since that is the input to COBRA
        k_row = k_samples,
        k_col = k_features,
        approximate = approx_neighbors,
        hnsw_args = hnsw_args
      )

      row_fusion <- fusion_wts$W_row$fusion
      col_fusion <- fusion_wts$W_col$fusion
      if(min(fusion_wts$W_row$values) <= 0 || min(fusion_wts$W_col$values) <= 0) {
        warning(paste(
          "U has diverged, try increasing k_features, k_samples",
          "arguments to encourage fusion terms.",
          "Returning current step"))
        rss_vals[t] <- sum(sweep(X - U, 2, sqrt(w ^ 2 + lambda * w), "*") ^ 2)
        objective_vals[t] <-
          gamma * (row_fusion + col_fusion) + rss_vals[t] / 2
        break
      }
    } else {
      row_fusion <- calc_fusion_term(fusion_wts$W_row, U, rows=TRUE)
      col_fusion <- calc_fusion_term(fusion_wts$W_col, U, rows=FALSE)
    }
    col_fusion_vals[t] <- col_fusion
    row_fusion_vals[t] <- row_fusion

    rss_vals[t] <- sum(sweep(X - U, 2, sqrt(w ^ 2 + lambda * w), "*") ^ 2)
    objective_vals[t] <- gamma * (row_fusion + col_fusion) + rss_vals[t] / 2

    if(which.min(objective_vals) == t) {
      best_obj_w <- w
      U_for_best_obj <- U
    }

    if(t > 1 && greedy_terminate && objective_vals[t] > objective_vals[t-1]) {
      break
    }

    cobra_result <-
      cobra(U_step,
            gamma = gamma,
            W_row = fusion_wts$W_row,
            W_col = fusion_wts$W_col,
            max_iter = tmax_cobra,
            max_iter_inner = 5 * tmax_cobra,
            tol = tol)

    U_prime <- cobra_result
    U_diff <- sum(abs(U_prime - U)) / sum(abs(U_prime))
    cobra_diffs[t] <- U_diff
    if (t > 1 && U_diff < tol) {
      U <- U_prime
      break
    }
    U <- U_prime

    col_sum_sq <- colSums((X - U) ^ 2)
    col_residuals[t, ] <- col_sum_sq

    lipschitz_fixed_U = sqrt(sum(col_sum_sq ^ 2))
    nu_w = 1/(1.1 * lipschitz_fixed_U)

    w_prime <- projection_onto_simplex(w - nu_w * (w + lambda / 2) * col_sum_sq)

    w_diff <- sum((w_prime - w) ^ 2)
    w_diffs[t] <- w_diff
    w <- w_prime
    w_path[t, ] <- w
    if (progress) {
      setTxtProgressBar(pb, t)
    }
  }

  w <- w_coordinate_descent(col_sum_sq, lambda)

  end <- Sys.time()

  if (progress) {
    close(pb)
  }

  to_return <- list(
    U = U,
    w = w,
    lambda = lambda,
    gamma = gamma,
    k_features = k_features,
    k_samples = k_samples,
    recalculate_weights = recalculate_weights,
    cobra_diffs = cobra_diffs,
    w_diffs = w_diffs,
    w_path = w_path,
    objective = objective_vals,
    rss = rss_vals,
    row_fusion = row_fusion_vals,
    col_fusion = col_fusion_vals,
    time = as.numeric(end - start, units = "secs"),
    col_residuals = col_residuals,
    best_obj_w = best_obj_w,
    U_for_best_obj = U_for_best_obj
  )

  to_return
}
