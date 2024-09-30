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
                         tmax_inner = 50,
                         tmax_outer = 100,
                         tol = 1e-4,
                         recalculate_weights = FALSE,
                         approx = 0,
                         progress = FALSE,
                         fusion_weights = NA,
                         return_fit = FALSE) {
  n <- nrow(X)
  p <- ncol(X)

  valid_indices <- which(is.finite(X), arr.ind = TRUE)
  invalid_indices <- which(!is.finite(X), arr.ind = TRUE)

  missing_columns <- sort(unique(invalid_indices[, 2]))
  stopifnot("missing values must have length > 0"= nrow(invalid_indices) > 0)

  filled_vals = matrix(NA, nrow=tmax_outer, ncol=nrow(invalid_indices))
  w_path = matrix(NA, nrow=tmax_outer, ncol=p)
  U = X
  U[invalid_indices] = mean(X[valid_indices])
  w = rep(1/p, p)

  if(!is.list(fusion_weights)) {
    fusion_wts <- fast_gkn_weights(
      t(U),
      k_col = k_samples,
      k_row = k_features,
      phi = phi,
      approx = approx
    )
  }


  for(t in 1:tmax_outer) {
    if (progress) {
      pb <- progress_bar$new(
        format = paste0(t, "/", tmax_outer, " imputing iter [:bar] :percent"),
        total = tmax_inner)
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
      tmax = tmax_inner,
      tol = tol
    )
    U <- bcbc_result$U
    w <- bcbc_result$w

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
                              tmax = NA,  # Set this to set both below
                              tmax_cobra = 100,
                              tmax_outer = 100,
                              tol = 1e-4,
                              recalculate_weights = FALSE,
                              approx = 0,
                              progress_bar = NA,
                              wts = NULL) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  U <- X
  w <- rep(1/p, p)

  if (!is.na(tmax)) {
    tmax_cobra <- tmax
    tmax_outer <- tmax
  }

  start <- Sys.time()
  for (t in 1:tmax_outer) {
    U_step <- U - sweep(X - U, 2, w ^ 2 + lambda * w, "*")
    # U_step <- U - sweep(X - U, 2, w, "*")

    U_step[missing_indices] <- U[missing_indices] -
      (X[missing_indices] - U[missing_indices]) *
      (w[missing_indices[, 2]] - w_tilde[missing_indices[, 2]])^2

    if (recalculate_weights || (t == 1 && !is.list(wts))) {
      wts <- fast_gkn_weights(
        t(U_step),
        k_col = k_samples,
        k_row = k_features,
        phi = phi,
        approx = approx
      )

      if(min(wts$w_row) <= 0 || min(wts$w_col) <= 0) {
        warning(paste(
          "U has diverged, try increasing k_samples, k_features or decreasing phi",
          "arguments to encourage fusion terms.",
          "Returning current step"))
        break
      }
    }

    # Cobra has rows as features, columns as samples
    cobra_result <-
      cobra(t(U_step),
            wts$E_row,
            wts$E_col,
            wts$w_row,
            wts$w_col,
            gamma = gamma,  # Need to double/half gamma?
            max_iter = tmax_cobra,
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
      lambda = lambda,
      time = as.numeric(end - start)
    )
  )
}

