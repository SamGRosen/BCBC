#' Perform BCBC with imputation of missing values
#'
#' @inheritParams BCBC
#' @param tmax_hierarchy vector of 3 values
#'  1. number of BCBC fits to run
#'  2. number of PALM iterations per BCBC fit
#'  3. number of COBRA iterations per PALM iteration
#' @param return_fits return all intermediate BCBC fits
#'
#' @return list with info
#'   1. `bcbc_fits` list of `BCBC` fits if `return_fits` is true
#'   2. `w` final fitted weight vector
#'   3. `w_path` fitted weight vectors from each BCBC fit
#'   4. `valid_rss` unweighted RSS of non-missing values
#'   5. `valid_weighted_rss` weighted RSS of non-missing values
#'   6. `filled_vals` imputed values from each BCBC fit
#'   7. `invalid_indices` indices of missing values
#'   8. `time` to run
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = FALSE)
#' missing_X <- checker$X
#' mask <- which(
#'   matrix(
#'     runif(nrow(checker$X) * ncol(checker$X)),
#'     nrow = nrow(checker$X)) < 0.2,
#'   arr.ind = TRUE)
#' heldout_vals <- checker$X[mask]
#' missing_X[mask] <- NA
#' missing_fits <- BCBC_missing(missing_X,
#'                              lambda = 0,
#'                              gamma = 200,
#'                              tmax_hierarchy = c(5, 50, 100),
#'                              recalculate_weights = TRUE,
#'                              return_fits = TRUE)
#' image(t(missing_X))
#' image(t(missing_fits$bcbc_fits[[5]]$U))
#' for(i in 1:5) {
#'   print(sum((missing_fits$filled_vals[i, ] - heldout_vals)^2))
#' }
BCBC_missing <- function(X,
                         lambda,
                         k_samples = 4,
                         k_features = 4,
                         gamma = 10,
                         phi = 1,
                         tmax_hierarchy = c(10, 25, 50),
                         tol = 1e-4,
                         recalculate_weights = FALSE,
                         approx_neighbors = FALSE,
                         hnsw_args = list(),
                         progress = FALSE,
                         fusion_wts = NA,
                         return_fits = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  valid_indices <- which(is.finite(X), arr.ind = TRUE)
  invalid_indices <- which(!is.finite(X), arr.ind = TRUE)

  missing_columns <- sort(unique(invalid_indices[, 2]))
  stopifnot("missing values must have length > 0"= nrow(invalid_indices) > 0)
  stopifnot("tmax_hierarchy must have length 3"= length(tmax_hierarchy) == 3)

  filled_vals = matrix(NA, nrow=tmax_hierarchy[1], ncol=nrow(invalid_indices))
  valid_rss <- rep(NA, tmax_hierarchy[1])
  valid_weighted_rss <- rep(NA, tmax_hierarchy[1])
  w_path = matrix(NA, nrow=tmax_hierarchy[1], ncol=p)
  bcbc_fits = list()

  U = X
  U[invalid_indices] = mean(X[valid_indices])
  w = rep(1/p, p)

  start <- Sys.time()

  for(t in 1:tmax_hierarchy[1]) {
    if (progress) {
      print(t)
    }

    M <- X
    M[invalid_indices] <- U[invalid_indices]
    bcbc_result <- BCBC(
      M,
      lambda = lambda,
      gamma = gamma,
      recalculate_weights = recalculate_weights,
      fusion_wts = fusion_wts,
      tmax_outer = tmax_hierarchy[2],
      tmax_cobra = tmax_hierarchy[3],
      tol = tol,
      hnsw_args = hnsw_args,
      approx_neighbors = approx_neighbors,
      init_U = M
    )

    if(return_fits) {
      bcbc_fits[[t]] <- bcbc_result
    }

    U <- bcbc_result$U
    w <- bcbc_result$w

    valid_rss[t] <- sum((X[valid_indices] - U[valid_indices]) ^ 2)
    swept_weighted_residuals <- sweep(X - U, 2, sqrt(w ^ 2 + lambda * w), "*") ^ 2
    valid_weighted_rss[t] <- sum(swept_weighted_residuals[valid_indices])
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

  end <- Sys.time()

  return(
    list(
      bcbc_fits = bcbc_fits,
      w = w,
      w_path = w_path,
      valid_rss = valid_rss,
      valid_weighted_rss = valid_weighted_rss,
      filled_vals = filled_vals,
      invalid_indices = invalid_indices,
      time = as.numeric(end - start, units = "secs")
    )
  )
}
