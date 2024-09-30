#' Generate simulated data resembling a checkerboard
#'
#' @param n number of rows
#' @param p number of columns
#' @param num_row_clusters number of row clusters
#' @param num_col_clusters number of column clusters
#' @param p_extra number of extra noisy columns
#' @param noise standard deviation of normal noise added to checkerboard
#' @param prob_empty between 0 and 1, probability a bicluster is 0
#' @param shuffle randomly reorder rows and columns
#' @param scale the resulting matrix
#'
#' @return
#' @export
#'
#' @examples
gen_checkerboard <- function(n,
                             p,
                             num_row_clusters,
                             num_col_clusters,
                             p_extra = 0,
                             noise = 1,
                             cluster_spread = 5,
                             prob_empty = 0.1,
                             shuffle = TRUE,
                             scale = TRUE) {
  row_partition <- sort(sample(1:num_row_clusters, n, replace = TRUE))
  col_partition <-
    sort(sample(1:num_col_clusters, p, replace = TRUE))

  mu_kr <- matrix(
    runif(num_row_clusters * num_col_clusters, -cluster_spread, cluster_spread),
    nrow = num_row_clusters,
    ncol = num_col_clusters
  )

  zeroed <- matrix(
    rbinom(num_row_clusters * num_col_clusters, 1, 1-prob_empty),
    nrow = num_row_clusters,
    ncol = num_col_clusters
  )
  mu_kr <- mu_kr * zeroed

  data_mat <- matrix(NA, nrow = n, ncol = (p + p_extra))
  centers <- matrix(NA, nrow = n, ncol = (p + p_extra))
  for (i in 1:n) {
    for (j in 1:p) {
      data_mat[i, j] <-
        rnorm(1, mean = mu_kr[row_partition[i], col_partition[j]],
              sd = noise)
      centers[i, j] <- mu_kr[row_partition[i], col_partition[j]]
    }
  }

  true_features = rep(1, p)
  if (p_extra > 0) {
    data_mat[, (p + 1):(p + p_extra)] <-
      rnorm(n * p_extra, mean = 0, sd = noise)
    centers[, (p + 1):(p + p_extra)] <- 0
    col_partition <- c(col_partition, rep(num_col_clusters + 1, p_extra))
    true_features <- c(true_features, rep(0, p_extra))
  }

  if (shuffle) {
    col_reorder <- sample(p + p_extra)
    row_reorder <- sample(n)
    shuffled_col <- data_mat[, col_reorder]
    shuffled_center <- centers[, col_reorder]
    shuffled_row <- shuffled_col[row_reorder,]
    shuffled_center <- centers[row_reorder,]
    row_partition <- row_partition[row_reorder]
    col_partition <- col_partition[col_reorder]
    true_features <- true_features[col_reorder]
    if(scale) {
      shuffled_row <- scale(shuffled_row)
      shuffled_center <- scale(shuffled_center,
                               center=attr(shuffled_row, "scaled:center"),
                               scale=attr(shuffled_row, "scaled:scale")
      )
    }
    return(list(X = shuffled_row, centers = shuffled_center,
                row_partition=row_partition, col_partition=col_partition,
                zeroed=zeroed, true_features=true_features))
  }
  if(scale) {
    data_mat <- scale(data_mat)
    centers <- scale(
      centers,
      center = attr(data_mat, "scaled:center"),
      scale = attr(data_mat, "scaled:scale")
    )
  }
  return(list(X = data_mat, centers = centers,
              row_partition=row_partition, col_partition=col_partition,
              zeroed=zeroed, true_features))
}
