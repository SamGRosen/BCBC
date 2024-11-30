library(tidyverse)
library(BCBC)
library(clusterSim)

generate_checkers <- function(seed=2024) {
  set.seed(seed)
  all_checkers <- list()

  trials <- 16
  noise_levels <- c(8)
  extra_dims <- c(0, 100, 200, 300, 400, 500, 600, 700, 800, 900)

  i <- 1
  for (trial in 1:trials) {
    for (noise_level in noise_levels) {
      for (extra_dim in extra_dims) {
        all_checkers[[i]] <- gen_checkerboard(
          200,
          200,
          5,
          5,
          noise = noise_level,
          cluster_spread = 10,
          p_extra = extra_dim,
          prob_empty = 0
        )
        all_checkers[[i]]$noise_level = noise_level
        all_checkers[[i]]$extra_dim = extra_dim
        all_checkers[[i]]$trial = trial
        i <- i + 1
      }
    }
  }

  all_checkers
}

generate_checkers_SNR <- function(seed=2024) {
  set.seed(seed)
  all_checkers <- list()

  trials <- 16
  noise_levels <- c(5, 7.5, 10, 12.5, 15)
  dimensions <- c(50, 75, 100, 125, 150)
  i <- 1
  for (trial in 1:trials) {
    for (noise_level in noise_levels) {
      for (dimension in dimensions) {
        all_checkers[[i]] <- gen_checkerboard(
          dimension,
          dimension,
          5,
          5,
          noise = noise_level,
          cluster_spread = 10,
          p_extra = 300,
          prob_empty = 0
        )
        all_checkers[[i]]$noise_level = noise_level
        all_checkers[[i]]$extra_dim = 300
        all_checkers[[i]]$trial = trial
        i <- i + 1
      }
    }
  }

  all_checkers
}

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
    row_result <- safe_compare(row_groups, checker$row_partition, metric)
    col_result <- safe_compare(col_groups, checker$col_partition, metric)
    col_result_true <- safe_compare(col_groups[true_feature_indices], checker$col_partition[true_feature_indices], metric)

    to_return[[paste(metric, "row", sep = "_")]] <- row_result
    to_return[[paste(metric, "col", sep = "_")]] <- col_result
    to_return[[paste(metric, "col_true", sep = "_")]] <- col_result_true
  }

  if(is.null(fitted_mat)) {
    to_return$true_err <- NA
    to_return$row_silhouette <- NA
    to_return$col_silhouette <- NA
    to_return$row_db <- NA
    to_return$col_db <- NA
    to_return$row_ch <- NA
    to_return$col_ch <- NA
  } else {
    to_return$true_err <- mean((checker$centers - fitted_mat)^2)
    to_return$row_silhouette <- index.S(dist(fitted_mat), checker$row_partition)
    to_return$col_silhouette <- index.S(dist(t(fitted_mat)), checker$col_partition)
    to_return$row_db <- index.DB(fitted_mat, checker$row_partition)$DB
    to_return$col_db <- index.DB(t(fitted_mat), checker$col_partition)$DB
    to_return$row_ch <- index.G1(fitted_mat, checker$row_partition)
    to_return$col_ch <- index.G1(t(fitted_mat), checker$col_partition)
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
  best_run <- bcbc_result$bcbc_cv$cv_data |> slice(which.min(BIC))
  best_index <- best_run$index
  list(
    row_groups = bcbc_result$bcbc_cv$row_clusters[[best_index]],
    col_groups = bcbc_result$bcbc_cv$col_clusters[[best_index]],
    fitted_mat = bcbc_result$bcbc_cv$all_runs[[best_index]]$U,
    fitted_feature_coefs = bcbc_result$bcbc_cv$all_runs[[best_index]]$w
  )
}

cobra_extractor <- function(cobra_result) {
  cobra_clusters = get_weighted_biclusters(t(cobra_result$U[[1]]),
                                           weights = NA,
                                           percent_noise = 0.5)
  list(
    col_groups = cobra_clusters$col_clusters,
    row_groups = cobra_clusters$row_clusters,
    fitted_mat = t(cobra_result$U[[1]])
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

