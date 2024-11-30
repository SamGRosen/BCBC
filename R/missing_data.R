#' Title
#'
#' @param X
#' @param lambda
#' @param k_samples
#' @param k_features
#' @param gamma
#' @param phi
#' @param tmax_inner
#' @param tmax_outer
#' @param tol
#' @param recalculate_weights
#' @param approx
#' @param progress
#'
#' @return
#' @import progress
#' @export
#'
#' @examples
BCBC_missing <- function(X,
                         lambda,
                         k_samples = 4,
                         k_features = 4,
                         gamma = 10,
                         phi = 0.05,
                         tmax_hierarchy = c(10, 25, 50),
                         tol = 1e-4,
                         recalculate_weights = FALSE,
                         approx_neighbors = FALSE,
                         hnsw_args = list(),
                         progress = FALSE,
                         fusion_weights = NA,
                         return_fit = FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  valid_indices <- which(is.finite(X), arr.ind = TRUE)
  invalid_indices <- which(!is.finite(X), arr.ind = TRUE)

  missing_columns <- sort(unique(invalid_indices[, 2]))
  stopifnot("missing values must have length > 0"= nrow(invalid_indices) > 0)
  stopifnot("tmax_hierarchy must have length 3"= length(tmax_hierarchy) == 3)

  filled_vals = matrix(NA, nrow=tmax_hierarchy[1], ncol=nrow(invalid_indices))
  row_fusions <- rep(NA, tmax_hierarchy[1])
  col_fusions <- rep(NA, tmax_hierarchy[1])
  valid_rss <- rep(NA, tmax_hierarchy[1])
  valid_weighted_rss <- rep(NA, tmax_hierarchy[1])
  w_path = matrix(NA, nrow=tmax_hierarchy[1], ncol=p)

  U = X
  U[invalid_indices] = mean(X[valid_indices])
  w = rep(1/p, p)

  if(!is.list(fusion_weights)) {
    fusion_wts <- fast_gkn_weights(
      t(U),
      k_col = k_samples,
      k_row = k_features,
      phi = phi,
      approximate = approx_neighbors,
      hnsw_args = hnsw_args
    )
  }


  for(t in 1:tmax_hierarchy[1]) {
    if (progress) {
      pb <- progress_bar$new(
        format = paste0(t, "/", tmax_hierarchy[1], " imputing iter [:bar] :percent"),
        total = tmax_hierarchy[2])
    } else {
      pb <- NULL
    }

    M <- X
    M[invalid_indices] <- U[invalid_indices]
    bcbc_result <- BCBC_missing_iter(
      M,
      w_tilde = w,
      missing_indices = invalid_indices,
      lambda = lambda,
      gamma = gamma,
      recalculate_weights = FALSE,
      wts = fusion_wts,
      progress_bar = pb,
      tmax_hierarchy = tmax_hierarchy[2:3],
      tol = tol
    )
    U <- bcbc_result$U
    w <- bcbc_result$w

    valid_rss[t] <- sum((X[valid_indices] - U[valid_indices]) ^ 2)
    swept_weighted_residuals <- sweep(X - U, 2, sqrt(w ^ 2 + lambda * w), "*") ^ 2
    valid_weighted_rss[t] <- sum(swept_weighted_residuals[valid_indices])
    row_fusions[t] <- bcbc_result$row_fusion
    col_fusions[t] <- bcbc_result$col_fusion
    filled_vals[t, ] <- U[invalid_indices]
    w_path[t, ] <- w

    if(t > 1 &&
       sum(abs(w_path[t - 1,] - w_path[t,])) < tol &&
       sum(abs(filled_vals[t - 1,] - filled_vals[t,])) / sum(abs(filled_vals[t, ])) < tol) {

      filled_vals <- filled_vals[1:t, ]
      w_path <- w_path[1:t, ]
      break
    }
  }

  return(
    list(
      U = ifelse(return_fit, U, NA),
      w = w,
      w_path = w_path,
      row_fusions = row_fusions,
      col_fusions = col_fusions,
      valid_rss = valid_rss,
      valid_weighted_rss = valid_weighted_rss,
      objectives = gamma * (row_fusions + col_fusions) + valid_weighted_rss / 2,
      filled_vals = filled_vals,
      invalid_indices = invalid_indices
    )
  )
}

BCBC_missing_iter <- function(X,
                              w_tilde,
                              missing_indices,
                              lambda,
                              k_samples = 4,
                              k_features = 4,
                              gamma = 10,
                              phi = 0.05,
                              tmax_hierarchy = c(25, 50),
                              tol = 1e-4,
                              recalculate_weights = FALSE,
                              approx_neigbors = FALSE,
                              hnsw_args = FALSE,
                              progress_bar = NA,
                              wts = NULL) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  U <- X
  w <- rep(1/p, p)

  stopifnot("tmax_hierarchy must have length 2"= length(tmax_hierarchy) == 2)

  start <- Sys.time()
  for (t in 1:tmax_hierarchy[1]) {
    U_step <- U - sweep(X - U, 2, w ^ 2 + lambda * w, "*")

    U_step[missing_indices] <- U[missing_indices] -
      (X[missing_indices] - U[missing_indices]) *
      (w[missing_indices[, 2]] - w_tilde[missing_indices[, 2]])^2

    if (recalculate_weights || (t == 1 && !is.list(wts))) {
      wts <- fast_gkn_weights(
        t(U_step),
        k_col = k_samples,
        k_row = k_features,
        phi = phi,
        approximate = approx_neighbors,
        hnsw_args = hnsw_args
      )
      row_fusion <- wts$row_fusion
      col_fusion <- wts$col_fusion

      if(min(wts$w_row) <= 0 || min(wts$w_col) <= 0) {
        warning(paste(
          "U has diverged, try increasing k_samples, k_features or decreasing phi",
          "arguments to encourage fusion terms.",
          "Returning current step"))
        break
      }
    } else {
      row_fusion <- calc_fusion_term(wts$unique_row_edges, wts$w_row, U, rows=FALSE)
      col_fusion <- calc_fusion_term(wts$unique_col_edges, wts$w_col, U, rows=TRUE)
    }

    # Cobra has rows as features, columns as samples
    cobra_result <-
      cobra(t(U_step),
            wts$E_row,
            wts$E_col,
            wts$w_row,
            wts$w_col,
            gamma = gamma,  # Need to double/half gamma?
            max_iter = tmax_hierarchy[2],
            tol=tol)

    U_prime <- t(cobra_result$U[[1]])
    U_diff <- sum(abs(U_prime - U)) / sum(abs(U_prime))
    if (U_diff < tol) {
      U <- U_prime
      break
    }
    U <- U_prime

    # w_deriv_as_mat <- matrix(1, nrow=n, ncol=p, byrow=TRUE)
    w_deriv_as_mat <- matrix(w + lambda / 2, nrow=n, ncol=p, byrow=TRUE)
    w_deriv_as_mat[missing_indices] <- w[missing_indices[, 2]] - w_tilde[missing_indices[, 2]]

    col_sum_sq <- colSums((X - U) ^ 2 * w_deriv_as_mat)
    lipschitz_fixed_U = sqrt(sum(col_sum_sq ^ 2))
    nu_w = 1/(1.1 * lipschitz_fixed_U)

    # w_prime <-
    #   projection_onto_simplex(w - nu_w / sqrt(p) * col_sum_sq)
    w_prime <- projection_onto_simplex(w - nu_w * col_sum_sq)
    w_diff <- sum((w_prime - w) ^ 2)
    w <- w_prime

    if (is.object(progress_bar) && t %% 10 == 0) {
      progress_bar$update(t / tmax_outer)
    }
  }

  end <- Sys.time()

  return(
    list(
      U = U,
      w = w,
      row_fusion = row_fusion,
      col_fusion = col_fusion,
      lambda = lambda,
      time = as.numeric(end - start)
    )
  )
}
