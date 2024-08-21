library(tidyverse)
library(BCBC)
library(WeightedCluster)

overlap_to_groups <- function(row_memberships, col_memberships) {
  row_membership_df <- data.frame(row_memberships)
  unique_row_groups <- row_membership_df |>
    group_by_all() |>
    summarise(COUNT = n(), .groups = 'drop') |>
    mutate(group = row_number())
  row_groups <- row_membership_df |>
    left_join(unique_row_groups, by = colnames(row_membership_df)) |>
    dplyr::select(group)

  col_membership_df <- data.frame(col_memberships)
  unique_col_groups <- col_membership_df |>
    group_by_all() |>
    summarise(COUNT = n(), .groups = 'drop') |>
    mutate(group = row_number())
  col_groups <- col_membership_df |>
    left_join(unique_col_groups, by = colnames(col_membership_df)) |>
    dplyr::select(group)

  list(
    row_groups = row_groups$group,
    col_groups = col_groups$group,
    fitted_mat = NULL
  )
}

metrics <- c("vi", "nmi", "split.join", "rand", "adjusted.rand")
evaluate_checker <- function(checker,
                             fitted_mat,
                             row_groups,
                             col_groups,
                             fitted_feature_coefs = NA) {
  to_return <- list()
  to_return$row_group_count <- length(unique(row_groups))
  to_return$col_group_count <- length(unique(col_groups))

  true_feature_indices <- which(checker$true_features > 0)
  untrue_feature_indices <- which(checker$true_features == 0)

  if(any(is.na(fitted_feature_coefs))) {
    to_return$feature_recall <- NA
    to_return$feature_precision <- NA
    to_return$feature_f1 <- NA
    to_return$feature_auc <- NA
  } else {
    TP <- sum(checker$true_features > 0 & fitted_feature_coefs > 0)
    FN <- sum(checker$true_features > 0 & fitted_feature_coefs == 0)
    FP <- sum(checker$true_features == 0 & fitted_feature_coefs > 0)
    to_return$feature_recall <- TP / (TP + FN)
    to_return$feature_precision <- TP / (TP + FP)
    if(to_return$feature_recall == 0 || to_return$feature_precision == 0) {
      to_return$feature_f1 <- 0
    } else {
      to_return$feature_f1 = 2 * TP / (2 * TP + FP + FN)
    }

    if(length(unique(checker$true_features)) == 2) {
      to_return$feature_auc <- ROCR::performance(
        ROCR::prediction(
          fitted_feature_coefs,
          checker$true_features), "auc")@y.values[[1]]
    } else {
      to_return$feature_auc <- NA
    }

  }

  if(length(untrue_feature_indices) > 0) {
    to_return$largest_group_in_untrue <- max(table(col_groups[untrue_feature_indices])) / length(untrue_feature_indices)
  } else {
    to_return$largest_group_in_untrue <- NA
  }
  for (metric in metrics) {
    if (length(row_groups) != length(checker$row_partition) || to_return$row_group_count == length(checker$row_partition)) {
      row_result <- -1
    } else {
      row_result <- igraph::compare(row_groups, checker$row_partition, method = metric)
    }
    if (length(col_groups) != length(checker$col_partition) || to_return$col_group_count == length(checker$col_partition)) {
      col_result <- -1
      col_result_true <- -1
    } else {
      col_result <- igraph::compare(col_groups, checker$col_partition, method = metric)
      col_result_true <- igraph::compare(col_groups[true_feature_indices],
                                         checker$col_partition[true_feature_indices],
                                         method = metric)
    }

    to_return[[paste(metric, "row", sep = "_")]] <- row_result
    to_return[[paste(metric, "col", sep = "_")]] <- col_result
    to_return[[paste(metric, "col_true", sep = "_")]] <- col_result_true
  }

  if(is.null(fitted_mat)) {
    to_return$true_err <- NA
  } else {
    to_return$true_err <- mean((checker$centers - fitted_mat)^2)
  }
  to_return
}

get_all_checker_results <- function(all_checkers, all_results, extractor) {
  all_info <- NULL
  for(run in 1:length(all_checkers)) {
    evaluated_info <- do.call(evaluate_checker,
                              c(list(checker = all_checkers[[run]]), extractor(all_results[[run]])))
    new_row <- as_tibble_row(
      c(evaluated_info,
        noise_level = all_checkers[[run]]$noise_level,
        trial = all_checkers[[run]]$trial,
        extra_dim = all_checkers[[run]]$extra_dim,
        true_row_num = length(unique(all_checkers[[run]]$row_partition)),
        true_col_num = length(unique(all_checkers[[run]]$col_partition)),
        n = length(all_checkers[[run]]$row_partition),
        p = length(all_checkers[[run]]$col_partition)
      )
    )
    if(is.null(all_info)) {
      all_info <- new_row
    } else {
      all_info <- all_info |> add_row(new_row)
    }
  }

  return(all_info)
}

bcbc_extractor <- function(bcbc_result) {
  best_run <- bcbc_result$cv_data |> slice(which.min(eBIC2))
  best_index <- best_run$index
  list(
    row_groups = bcbc_result$row_clusters[[best_index]],
    col_groups = bcbc_result$col_clusters[[best_index]],
    fitted_mat = bcbc_result$all_runs[[best_index]]$U,
    fitted_feature_coefs = bcbc_result$all_runs[[best_index]]$w
  )
}

bcbc_extractor2 <- function(bcbc_result) {
  swept <- sweep(bcbc_result$bcbc$U,
                 2,
                 bcbc_result$bcbc$w ^ 2 + bcbc_result$bcbc$lambda * bcbc_result$bcbc$w,
                 "*")

  row_threshold <- get_threshold_for_k_components(dist(swept),
                                                    5)

  row_clusters <- get_row_clusters(swept, row_threshold)
  col_clusters <- wcKMedoids(dist(t(bcbc_result$bcbc$U)), 6, weights = bcbc_result$bcbc$w)$clustering
  list(
    row_groups = row_clusters,
    col_groups = col_clusters,
    fitted_mat = bcbc_result$bcbc$U,
    fitted_feature_coefs = bcbc_result$bcbc$w
  )
}

cobra_extractor <- function(cobra_result) {
  best_result <- which.min(cobra_result$validation_error)
  list(
    col_groups = cobra_result$groups_row[[best_result]]$cluster,
    row_groups = cobra_result$groups_col[[best_result]]$cluster,
    fitted_mat = t(cobra_result$U[[best_result]])
  )
}

bcel_extractor <- function(bcel_result) {
  to_return <- overlap_to_groups(abs(bcel_result$U) > 0, abs(bcel_result$V) > 0)
  to_return$fitted_mat <- bcel_result$U %*% t(bcel_result$V)
  return(to_return)
}

biclust_extractor <- function(biclust_obj) {
  overlap_to_groups(biclust_obj@RowxNumber, t(biclust_obj@NumberxCol))
}

sparseBC_extractor <- function(sparseBC_run) {
  list(
    row_groups = sparseBC_run$Cs,
    col_groups = sparseBC_run$Ds,
    fitted_mat = sparseBC_run$mus
  )
}

scbiclust_extractor <- function(scbiclust_result) {
  stacked_rows <- do.call(rbind, scbiclust_result$which.x)
  stacked_cols <- do.call(rbind, scbiclust_result$which.y)
  overlap_to_groups(t(stacked_rows), t(stacked_cols))
}

dct_extractor <- identity

