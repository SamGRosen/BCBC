#' Dykstra-like proximal algorithm from Convex Biclustering by Chi et al.
#'
#' @param X data matrix
#' @param gamma regularization parameter
#' @param W_row `sparseweights` object output from `knn_graph` for row fusion
#' @param W_col `sparseweights` object output from `knn_graph` for col fusion
#' @param max_iter number of convex cluster optimizations to solve
#' @param max_iter_inner number of inner iterations for `CCMMR`
#' @param tol of convergence
#'
#' @return fitted bicluster matrix
#' @import CCMMR
#' @seealso [fast_gkn_weights()]
#' @seealso [knn_graph()]
#' @seealso [knn_graph_approx()]
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 0, shuffle = FALSE)
#' W_row <- knn_graph(checker$X, 5)
#' W_col <- knn_graph(t(checker$X), 5)
#' cobra_fit <- cobra(checker$X, gamma = 10000, W_row = W_row, W_col = W_col)
#'
#' plot_fit(matrix_fit_to_df(cobra_fit))
cobra <- function(X,
                  gamma,
                  W_row,
                  W_col,
                  max_iter = 1e2,
                  max_iter_inner = 1e2,
                  tol = 1e-4) {
  U = matrix(X, nrow = nrow(X), ncol = ncol(X))
  P = matrix(0, nrow = nrow(X), ncol = ncol(X))
  Q = matrix(0, nrow = ncol(X), ncol = nrow(X))

  diff = 10 * tol
  iter = 1
  while (diff > tol && iter < max_iter) {
    Y = convex_clusterpath(
      t(U) + t(P),
      W_col,
      gamma,
      center = F,
      scale = F,
      max_iter_conv = max_iter_inner,
      eps_conv = tol
    )$coordinates

    P_prime = U + P - t(Y)

    U_prime = convex_clusterpath(
      t(Y) + t(Q),
      W_row,
      gamma,
      center = F,
      scale = F,
      max_iter_conv = max_iter_inner,
      eps_conv = tol
    )$coordinate
    Q_prime = Y + Q - t(U_prime)

    diff = sqrt(sum((U_prime - t(Y)) ^ 2)) / ncol(X) / nrow(X)

    P = P_prime
    U = U_prime
    Q = Q_prime
    iter = iter + 1
  }
  U
}


#' Utility method to use COBRA internals and calculate weights
#'
#' @param k_samples number of nearest neighbors for sample affinity graph
#' @param k_features number of nearest neighbors for feature affinity graph
#' @inherit cobra params
#' @inherit cobra return
#' @inherit cobra seealso
#' @return fitted bicluster matrix
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = FALSE)
#' cobra_knn(checker$X, 200, 4, 4)
cobra_knn <- function(X,
                      gamma,
                      k_samples = 2,
                      k_features = 2,
                      max_iter = 250,
                      max_iter_inner = 250,
                      tol = 1e-4) {
  wts <- fast_gkn_weights(
    X,
    k_row = k_samples,
    k_col = k_features,
    approx = FALSE
  )
  cobra(
    X,
    gamma = gamma,
    W_row = wts$W_row,
    W_col = wts$W_col,
    max_iter = max_iter,
    max_iter_inner = max_iter_inner,
    tol = tol
  )
}


#' Algorithm 2 of Convex Biclustering by Chi et al.
#'
#' @inheritParams cobra
#' @param tmax_hierarchy vector of 3 values
#'  1. number of COBRA fits
#'  2. number of number of Convex Clustering fits per COBRA fit
#'  3. number of iterations per Convex Clustering fit
#' @param progress show progress
#' @return list with info
#'   4. `valid_rss` unweighted RSS of non-missing values
#'   6. `filled_vals` imputed values from each COBRA fit
#'   7. `invalid_indices` indices of missing values
#'   8. `time` to run
#' @export
COBRA_missing <- function(X,
                          W_row,
                          W_col,
                          gamma = 10,
                          tmax_hierarchy = c(50, 10, 1000),
                          tol = 1e-4,
                          progress = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  valid_indices <- which(is.finite(X), arr.ind = TRUE)
  invalid_indices <- which(!is.finite(X), arr.ind = TRUE)

  missing_columns <- sort(unique(invalid_indices[, 2]))
  stopifnot("missing values must have length > 0"= nrow(invalid_indices) > 0)
  stopifnot("tmax_hierarchy must have length 3"= length(tmax_hierarchy) == 3)

  U = X
  U[invalid_indices] = mean(X[valid_indices])
  valid_rss <- rep(NA, tmax_hierarchy[1])
  filled_vals = matrix(NA, nrow=tmax_hierarchy[1], ncol=nrow(invalid_indices))

  start <- Sys.time()
  for(t in 1:tmax_hierarchy[1]) {
    if (progress) {
      print(t)
    }

    M <- X
    M[invalid_indices] <- U[invalid_indices]
    cobra_result <- cobra(
      M,
      gamma = gamma,
      W_row = W_row,
      W_col = W_col,
      max_iter = tmax_hierarchy[2],
      max_iter_inner = tmax_hierarchy[3],
      tol = tol
    )

    U <- cobra_result

    valid_rss[t] <- sum((X[valid_indices] - U[valid_indices]) ^ 2)

    filled_vals[t, ] <- U[invalid_indices]

    if(t > 1 &&
       sum(abs(filled_vals[t - 1,] - filled_vals[t,])) / sum(abs(filled_vals[t, ])) < tol) {

      filled_vals <- filled_vals[1:t, ]
      break
    }
  }

  end <- Sys.time()

  return(
    list(
      valid_rss = valid_rss,
      filled_vals = filled_vals,
      invalid_indices = invalid_indices,
      time = as.numeric(end - start, units = "secs")
    )
  )
}
