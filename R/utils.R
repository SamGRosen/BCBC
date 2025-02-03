#' Title
#'
#' @param estimated_clusters
#' @param true_clusters
#' @param method
#'
#' @return
#' @export
#' @importFrom igraph compare
#'
#' @examples
safe_compare <- function(estimated_clusters, true_clusters, method) {
  if (length(estimated_clusters) != length(true_clusters) ||
      length(unique(estimated_clusters)) == length(true_clusters)) {
    return(-1)
  } else {
    return(igraph::compare(estimated_clusters, true_clusters, method = method))
  }
}


#' Title
#'
#' @param matrix
#' @param labels
#' @param labels2
#' @param weights
#' @param lambda
#' @param filter_weight
#' @param cluster_with_weights
#'
#' @return
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom dplyr inner_join
#' @importFrom dplyr distinct
#' @importFrom tidyr pivot_longer
#' @export
#'
#' @examples
matrix_fit_to_df <- function(matrix,
                             labels=NULL,
                             labels2=NULL,
                             weights=c(),
                             lambda=0,
                             filter_weight=0,
                             cluster_with_weights=TRUE) {
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
  h_col <- hclust(col_dist)
  df <- df |>
    inner_join(tibble(
      row_num = h_row$order,
      order_row = 1:length(h_row$order)
    ), by = join_by(row_num)) |>
    inner_join(tibble(
      col_num = h_col$order,
      order_col = 1:length(h_col$order)
    ), by = join_by(col_num))

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


#' Title
#'
#' @param bcbc_result
#' @param cluster_with_weights
#' @param filter_weight
#' @param labels
#' @param labels2
#'
#' @return
#' @export
#'
#' @examples
bcbc_result_to_df <- function(bcbc_result,
                              cluster_with_weights = TRUE,
                              filter_weight = -Inf,
                              labels=NULL,
                              labels2=NULL) {
  # filter weight will be as if columns with weight below filter weight did not exist
  matrix_fit_to_df(
    bcbc_result$U,
    labels = labels,
    labels2 = labels2,
    weights = bcbc_result$w,
    lambda = bcbc_result$lambda,
    filter_weight = filter_weight,
    cluster_with_weights = cluster_with_weights
  )
}


#' Title
#'
#' @param plot_df
#' @param title
#' @param xlabel
#' @param ylabel
#' @param fill_attr
#' @param alpha_weight
#' @param bin_scale
#'
#' @return
#' @import ggplot2
#' @export
#'
#' @examples
plot_fit <- function(plot_df,
                     title = "",
                     xlabel = "",
                     ylabel = "",
                     fill_attr = "value",
                     alpha_weight = FALSE,
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


#' Title
#'
#' @param w_path
#' @param reorder
#'
#' @return
#' @export
#'
#' @examples
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
