#' Orthogonal projection of a vector onto the simplex
#'
#' From https://github.com/kharchenkolab/vrnmf/tree/main
#' @param unproj
#' @param bound
#'
#' @return
#' @export
#'
#' @examples
projection_onto_simplex <- function(unproj, bound=1) {
  q <- sort(unproj, decreasing = TRUE, method = "quick")
  qcum <- cumsum(q)
  mu <- (qcum - bound) / 1:length(qcum)
  cond1 <- (mu[-length(mu)] - q[-1]) > 0
  if (max(cond1) == 0) {
    ind <- length(mu)
  } else{
    ind <- which.max(cond1)
  }
  return(pmax(0, unproj - mu[ind]))
}


#' Title
#'
#' @param edges
#' @param weights
#' @param U
#' @param rows
#'
#' @return
#' @export
#'
#' @examples
calc_fusion_term <- function(edges, weights, U, rows=TRUE) {
  from_node <- edges$from_node
  to_node <- edges$to_node

  all_distances <- rep(NA, length(from_node))
  for(i in seq_along(from_node)) {
    if(rows) {
      all_distances[i] = sqrt(sum((U[from_node[i], ] - U[to_node[i], ])^2))
    } else {
      all_distances[i] = sqrt(sum((U[, from_node[i]] - U[, to_node[i]])^2))
    }
  }
  return(sum(weights * all_distances))
}


#' Title
#'
#' @param col_sum_sq Sum of squared errors for columns
#' @param lambda
#'
#' @return
#' @export
#'
#' @examples
w_coordinate_descent <- function(col_sum_sq, lambda) {
  f = function(alpha) {
    s = sum(pmax(alpha / col_sum_sq - lambda, 0))
    return(s / 2 - 1)
  }
  alpha = uniroot(f, c(0.0, 1000.0), tol = 1e-10)$root

  pmax(alpha / col_sum_sq - lambda, 0) / 2
}


#' Title
#'
#' @param X
#' @param lambda
#' @param k_features
#' @param k_samples
#' @param gamma
#' @param phi
#' @param tmax
#' @param tmax_cobra
#' @param tmax_outer
#' @param tol
#' @param recalculate_weights
#' @param greedy_terminate
#' @param approx_neighbors
#' @param hnsw_args
#' @param progress
#' @param fusion_wts
#'
#' @return
#' @export
#'
#' @examples
BCBC <- function(X,
                 lambda,
                 k_features = 4,
                 k_samples = 4,
                 gamma = 10,
                 phi = 1,
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
      t(X),
      k_row = k_features,
      k_col = k_samples,
      phi = phi,
      approximate = approx_neighbors,
      hnsw_args = hnsw_args
    )
  }

  for (t in 1:tmax_outer) {
    U_step <- U - sweep(U - X, 2, w ^ 2 + lambda * w, "*")

    if (recalculate_weights && t > 1) {
      fusion_wts <- fast_gkn_weights(
        t(U_step),  # After first step, use U_step here since that is the input to COBRA
        k_row = k_features,
        k_col = k_samples,
        phi = phi,
        approximate = approx_neighbors,
        hnsw_args = hnsw_args
      )

      row_fusion <- fusion_wts$row_fusion
      col_fusion <- fusion_wts$col_fusion
      if(min(fusion_wts$w_row) <= 0 || min(fusion_wts$w_col) <= 0) {
        warning(paste(
          "U has diverged, try increasing k_features, k_samples or decreasing phi",
          "arguments to encourage fusion terms.",
          "Returning current step"))
        rss_vals[t] <- sum(sweep(X - U, 2, sqrt(w ^ 2 + lambda * w), "*") ^ 2)
        objective_vals[t] <-
          gamma * (row_fusion + col_fusion) + rss_vals[t] / 2
        break
      }
    } else {
      row_fusion <- calc_fusion_term(fusion_wts$unique_row_edges, fusion_wts$w_row, U, rows=FALSE)
      col_fusion <- calc_fusion_term(fusion_wts$unique_col_edges, fusion_wts$w_col, U, rows=TRUE)
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

    # Cobra has rows as features, columns as samples
    cobra_result <-
      cobra(t(U_step),
            fusion_wts$E_row,
            fusion_wts$E_col,
            fusion_wts$w_row,
            fusion_wts$w_col,
            gamma = gamma,
            max_iter = tmax_cobra,
            tol=tol)

    U_prime <- t(cobra_result$U[[1]])
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
    # w_prime <- projection_onto_simplex(w - nu_w / sqrt(p) * col_sum_sq)

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
    phi = phi,
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
