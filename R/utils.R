#' Method to easily use `igraph::compare`
#'
#' @param estimated_clusters community memberships
#' @param true_clusters community memberships
#' @param method passed to `igraph::compare`
#'
#' @return metric in community memberships
#' @export
#' @importFrom igraph compare
#' @seealso [igraph::compare()]
safe_compare <- function(estimated_clusters, true_clusters, method) {
  if (length(estimated_clusters) != length(true_clusters) ||
      length(unique(estimated_clusters)) == length(true_clusters)) {
    return(-1)
  } else {
    return(igraph::compare(estimated_clusters, true_clusters, method = method))
  }
}


#' Convert a matrix into a dataframe for post-processing or plotting.
#'
#' @param matrix to convert
#' @param labels to match to each row
#' @param labels2 secondary set to match to each row
#' @param weights to measure importance of features
#' @param lambda hyperparameter used to create fit
#' @param filter_weight remove features with weight less than or equal to this
#' @param cluster_with_weights perform clustering with feature weights
#' @param method passed to `hclust` to determine row and column partitioning

#' @return data frame with column names
#'   1. `row_num` original row number of matrix entry
#'   1. `col_num` original column number of matrix entry
#'   1. `value` original value of matrix entry
#'   1. `order_row` a permutated row number for plotting similar rows together
#'   1. `order_col` a permutated column number for plotting similar columns together
#'   1. `weight` weight of original column if `weights` parameter provided
#'   1. `label` label of original row if `labels` parameter provided
#'   1. `label2` label of original row if `labels2` parameter provided
#'
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom dplyr inner_join
#' @importFrom dplyr distinct
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' ggplot(matrix_fit_to_df(checker$X), aes(x = order_col, y = order_row)) +
#'   geom_raster(aes(fill = value))
matrix_fit_to_df <- function(matrix,
                             labels=NULL,
                             labels2=NULL,
                             weights=c(),
                             lambda=0,
                             filter_weight=0,
                             cluster_with_weights=TRUE,
                             method = "ward.D") {
  if(length(weights) > 0) {
    valid_weights <- which(weights > filter_weight)
    weights <- weights[valid_weights]
    matrix <- matrix[, valid_weights]
  }
  df <- as.data.frame(matrix) |>
    mutate(row_num = row_number()) |>
    pivot_longer(!row_num, names_to = "col_num") |>
    mutate(col_num = as.numeric(sub("V", "", col_num)))

  if(length(weights > 0) && cluster_with_weights) {
    w_vals = sqrt(weights^2 + lambda * weights)
    matrix <- sweep(matrix, 2, w_vals, "*")
  }

  row_dist <- dist(matrix)
  col_dist <- dist(t(matrix))

  h_row <- hclust(row_dist)
  h_col <- hclust(col_dist, method = method)
  df <- df |>
    inner_join(tibble(
      row_num = h_row$order,
      order_row = 1:length(h_row$order)
    ), by = join_by(row_num)) |>
    inner_join(tibble(
      col_num = h_col$order,
      order_col = 1:length(h_col$order)
    ), by = join_by(col_num))

  if (length(weights) > 0) {
    df <- df |>
      inner_join(tibble(col_num = 1:length(valid_weights),
                        weight = weights[valid_weights]),
                 by = join_by(col_num))
  }

  if (is.null(labels)) {
    return(df)
  } else {
    df <- df |>
      inner_join(tibble(row_num = 1:length(h_row$order), label = labels),
                 by = join_by(row_num))
  }

  if (is.null(labels2)) {
    return(df)
  } else {
    df <- df |>
      inner_join(tibble(row_num = 1:length(h_row$order), label2 = labels2,
                        by = join_by(row_num)))
  }
  return(df)
}


#' Convert a BCBC result from `BCBC`.
#'
#' @param bcbc_result from `BCBC`
#' @inheritParams matrix_fit_to_df
#' @inherit matrix_fit_to_df return
#'
#' @export
#' @seealso [matrix_fit_to_df()]
#' @seealso [BCBC::BCBC()]
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.01, gamma = 100)
#' as_df <- bcbc_result_to_df(bcbc_fit)
#' ggplot(as_df, aes(x = order_col, y = order_row)) +
#'   geom_raster(aes(fill = value, alpha = weight))
bcbc_result_to_df <- function(bcbc_result,
                              cluster_with_weights = TRUE,
                              filter_weight = -Inf,
                              labels = NULL,
                              labels2 = NULL,
                              method = "ward.D") {
  # filter weight will be as if columns with weight below filter weight did not exist
  matrix_fit_to_df(
    bcbc_result$U,
    labels = labels,
    labels2 = labels2,
    weights = bcbc_result$w,
    lambda = bcbc_result$lambda,
    filter_weight = filter_weight,
    cluster_with_weights = cluster_with_weights,
    method = method
  )
}


#' Conveinence function for plotting matrix fits from `matrix_fit_to_df` and
#'   `bcbc_result_to_df`.
#'
#' @param plot_df return value of `matrix_fit_to_df` or `bcbc_result_to_df`.
#' @param title to display
#' @param xlabel x-axis label to display
#' @param ylabel y-axis label to display
#' @param fill_attr color the entries based off this column from `plot_df`
#' @param alpha_weight match opacity of columns with `weight` column from `plot_df`
#' @param weight_annotation add a line plot annotating weight magnitude of columns
#' @param bin_scale use a discrete but ordered color scale with this number of bins
#'
#' @return ggplot object
#' @import ggplot2
#' @export
#' @seealso [matrix_fit_to_df()]
#' @seealso [bcbc_result_to_df()]
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.01, gamma = 100)
#' as_df <- bcbc_result_to_df(bcbc_fit)
#' plot_fit(as_df, title="Title", xlabel="columns", ylabel="rows", bin_scale = 5)
#' plot_fit(as_df, title="Title", weight_annotation = TRUE)
plot_fit <- function(plot_df,
                     title = "",
                     xlabel = "",
                     ylabel = "",
                     fill_attr = "value",
                     alpha_weight = FALSE,
                     weight_annotation = FALSE,
                     bin_scale = FALSE) {
  if("label" %in% colnames(plot_df)) {
    label_annotation <- annotate(
      geom = "rect",
      xmin = -max(plot_df$order_col) / 25,
      xmax = 0,
      ymin = plot_df$order_row - 0.5,
      ymax = plot_df$order_row + 0.5,
      fill = case_when(
        plot_df$label == 0 ~ "#11DDDD",
        plot_df$label == 1 ~ "green",
        plot_df$label == 2 ~ "yellow",
        plot_df$label == 3 ~ "purple",
        plot_df$label == 4 ~ "grey",
        plot_df$label == 5 ~ "orange",
        plot_df$label == 6 ~ "blue",
        plot_df$label == 7 ~ "red",
        plot_df$label == 8 ~ "magenta",
        plot_df$label == 9 ~ "pink",
        TRUE ~ "black"
      )
    )
  } else {
    label_annotation <- geom_blank()
  }

  if(weight_annotation) {
    weight_annotation <- annotate(
      geom = "line",
      x = plot_df$order_col,
      y = (plot_df$weight - min(plot_df$weight)) * 10000
    )
  } else {
    weight_annotation <- geom_blank()
  }

  if (alpha_weight) {
    raster <- geom_raster(aes(fill = .data[[fill_attr]], alpha = weight))
  } else {
    raster <- geom_raster(aes(fill = .data[[fill_attr]]))
  }

  if (!bin_scale) {
    plot_scale = scale_fill_gradient2()
  } else {
    plot_scale = scale_fill_steps2(n.breaks = bin_scale)
  }

  ggplot(plot_df, aes(x = order_col, y = order_row)) +
    raster +
    label_annotation +
    weight_annotation +
    plot_scale +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      # axis.title.x = element_blank(),
      # axis.title.y = element_blank(),
      plot.title = element_text(size=30)
    ) +
    guides(fill = "none") +
    xlab(xlabel) +
    ylab(ylabel) +
    ggtitle(title)
}


#' Debug BCBC fits by examining the weights during optimization iteration
#'
#' @param w_path from `$w_path` of a `BCBC` fit value
#' @param reorder the columns based on weight values
#' @param extra_data matrix with equal dimensions as weight path from BCBC fit
#'  to plot instead
#'
#' @return ggplot_object
#' @export
#'
#' @examples
#' checker <- gen_checkerboard(100, 150, 5, 5, p_extra = 100, shuffle = TRUE)
#' bcbc_fit <- BCBC(checker$X, lambda = 0.05, gamma = 200)
#' plot_w_path(bcbc_fit$w_path)
#' plot_w_path(bcbc_fit$w_path, extra_data=bcbc_fit$col_residuals)
plot_w_path <- function(w_path, reorder=TRUE, extra_data=w_path) {
  if (any(is.na(w_path[, 1]))) {
    first_NA <- which.max(is.na(w_path[, 1]))
    w_path <- w_path[1:(first_NA - 1), ]
  }

  if(reorder) {
    w_path <- w_path[, order(colSums(w_path == 0), -colSums(w_path))]
  }
  as_df <- as.data.frame(extra_data) |>
    mutate(row_num = row_number()) |>
    pivot_longer(!row_num, names_to = "col_num") |>
    mutate(col_num = as.numeric(sub("V", "", col_num)))

  ggplot(as_df, aes(col_num, row_num)) +
    geom_raster(aes(fill = value)) +
    xlab("Index") +
    ylab("Iteration") +
    # scale_fill_gradient2(low = "#EEE", high = "blue") +
    scale_fill_gradient2(low="red", mid = "#EEE", high = "blue", midpoint = 1/ncol(w_path)) +
    annotate(
      "line",
      x = 1:ncol(w_path),
      y = colSums(w_path != 0),
      alpha = 0.25,
      linewidth = 1.1,
      color = "red"
    ) +
    annotate(
      "line",
      x = rowSums(w_path != 0),
      y = 1:nrow(w_path),
      orientation = "y",
      alpha = 0.25,
      linewidth = 1.1,
      color = "green"
    ) +
    theme_minimal()
}
