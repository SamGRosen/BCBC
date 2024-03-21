library(cvxbiclustr)
library(tidyverse)
library(dbscan)
library(future.apply)

solve_alpha_prox = function(D_l, nu, q, lambda) {
  D_nu <- D_l / (1 / nu + D_l)
  f <- function(alpha) {
    s <- sum(D_nu * pmax((alpha + q / nu) / D_l - lambda / 2, 0))
    return(s - 1)
  }
  a <- uniroot(f, c(-min(q / nu), 10000.0), tol = 1e-8)$root
  return(a)
}

create_edge_incidence_edges <- function(P, n) {
  nEdges <- nrow(P)
  E <- Matrix(
    data = 0,
    nrow = nEdges,
    ncol = n,
    sparse = TRUE
  )
  r <- 1:nEdges
  col <- P[, 1]
  E[(col - 1) * nEdges + r] <- 1
  col <- P[, 2]
  E[(col - 1) * nEdges + r] <- -1
  return(E)
}

fast_gkn_weights <- function(X, k_row, k_col, phi, approx = 0) {
  p <- ncol(X)
  n <- nrow(X)
  
  all_row_knn <- hnsw_knn(X, k_row + 1, byrow=FALSE)
  all_col_knn <- hnsw_knn(X, k_col + 1, byrow=FALSE)
  
  # all_row_knn <- FNN::get.knn(X, k_row, algorithm="cover_tree")
  # all_col_knn <- FNN::get.knn(t(X), k_col, algorithm="cover_tree")
  unique_row_edges <- data.frame(all_row_knn$idx) |>
    mutate(index = row_number()) |>
    pivot_longer(!index) |>
    mutate(
      neighbor_index = as.integer(str_remove(name, "X")),
      from_node = pmin(index, value),
      to_node = pmax(index, value)
    ) |>
    filter(from_node != to_node) |>
    distinct(from_node, to_node, .keep_all = TRUE)
  
  unique_col_edges <- data.frame(all_col_knn$idx) |>
    mutate(index = row_number()) |>
    pivot_longer(!index) |>
    mutate(
      neighbor_index = as.integer(str_remove(name, "X")),
      from_node = pmin(index, value),
      to_node = pmax(index, value)
    ) |>
    filter(from_node != to_node) |>
    distinct(from_node, to_node, .keep_all = TRUE)
  
  all_row_dist_sq <- all_row_knn$dist[cbind(unique_row_edges$from_node,
                                            unique_row_edges$neighbor_index)]^2
  all_col_dist_sq <- all_col_knn$dist[cbind(unique_col_edges$from_node,
                                            unique_col_edges$neighbor_index)]^2
  
  row_weights <- exp(-all_row_dist_sq * phi / p)
  row_weights <- row_weights / sum(row_weights) / sqrt(p)
  row_fusion <- sum(row_weights * all_row_dist_sq)
  col_weights <- exp(-all_col_dist_sq * phi / n)
  col_weights <- col_weights / sum(col_weights) / sqrt(n)
  col_fusion <- sum(col_weights * all_col_dist_sq)
  
  # Construct edge-incidence matrices
  E_row <-
    create_edge_incidence_edges(cbind(unique_row_edges$from_node, unique_row_edges$to_node),
                                n)
  E_col <-
    create_edge_incidence_edges(cbind(unique_col_edges$from_node, unique_col_edges$to_node),
                                p)
  
  return(list(
    w_row = row_weights,
    w_col = col_weights,
    E_row = E_row,
    E_col = E_col,
    row_fusion = row_fusion,
    col_fusion = col_fusion
  ))
}

rscobra <- function(X,
                    lambda,
                    k_row = 4,
                    k_col = 4,
                    nu = 1,
                    gamma = 10,
                    phi = 0.05,
                    tmax = NA,  # Set this to set all 3 below
                    tmax_cobra = 100,
                    tmax_biconvex = 100,
                    tmax_outer = 100,
                    tol = 1e-6,
                    recalculate_weights = TRUE,
                    approx = 0,
                    threshold = 0,
                    adaptive_nu = FALSE,  # Experimental, may not be helpful
                    progress = TRUE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  Q <- X
  q <- runif(p)
  q <- q / sum(q)
  
  if(!is.na(tmax)) {
    tmax_cobra <- tmax
    tmax_biconvex <- tmax
    tmax_outer <- tmax
  }
  
  if(adaptive_nu) {
    adaptive_steps <- accumulate(1:tmax_outer, function(x, unused) {1/2 + sqrt(1 + 4 * x^2)/2})
  }
  
  cobra_diffs <- rep(NA, tmax_outer)
  biconvex_diffs <- matrix(NA, nrow = tmax_outer, ncol = tmax_biconvex)
  w_path <- matrix(NA, nrow = tmax_outer, ncol = p)
  objective_vals <- rep(NA, tmax_outer)
  rss_vals <- rep(NA, tmax_outer)
  
  if(progress) {
    pb <- txtProgressBar(min = 0,
                         max = tmax_outer,
                         initial = 0)
  }
  
  start <- Sys.time()
  for (t in 1:tmax_outer) {
    if (recalculate_weights || t == 1) {
      # wts = gkn_weights(t(Q), k_row=k_row, k_col=k_col, phi=phi, return_connectivity=FALSE)
      wts <- fast_gkn_weights(
        t(Q),
        k_row = k_row,
        k_col = k_col,
        phi = phi,
        approx = approx
      )
      w_row <- wts$w_row
      w_col <- wts$w_col
      E_row <- wts$E_row
      E_col <- wts$E_col
      row_fusion <- wts$row_fusion
      col_fusion <- wts$col_fusion
    }
    
    rss_vals[t] <- sum(sweep(X - Q, 2, q ^ 2 + lambda * q, "*") ^ 2)
    objective_vals[t] <- gamma * (row_fusion + col_fusion) + rss_vals[t] / 2
    
    # Cobra has rows as features, columns as samples
    cobra_result <-
      cobra(t(Q),
            E_row,
            E_col,
            w_row,
            w_col,
            gamma = gamma,
            max_iter = tmax_cobra,
            tol=tol)
    
    Q_prime <- t(cobra_result$U[[1]])
    Q_diff <- sum(abs(Q_prime - Q)) / sum(abs(Q_prime))
    cobra_diffs[t] <- Q_diff
    if (Q_diff < tol) {
      Q <- Q_prime
      break
    }
    old_Q <- Q
    Q <- Q_prime
    
    for (t2 in 1:tmax_biconvex) {
      # Coordinate wise minima of biconvex optimization
      # Solve for Q with fixed w
      w_sq <- q ^ 2 + lambda * q
      X_weighted <- sweep(X, 2, w_sq / (1 / nu + w_sq), "*")
      Q_weighted <- sweep(Q, 2, 1 / (1 + nu * w_sq), "*")
      Q <- Q_weighted + X_weighted
      
      # Solve for w with fixed U
      col_sum_sq <- colSums((X - Q) ^ 2)
      alpha <- solve_alpha_prox(col_sum_sq, nu, q, lambda)
      new_q <- (col_sum_sq) / (col_sum_sq + 1 / nu) * pmax((alpha + q / nu) / col_sum_sq - lambda /
                                                             2, 0)
      
      q_diff <- sum((new_q - q) ^ 2)
      biconvex_diffs[t, t2] <- q_diff
      if (q_diff < tol) {
        q <- new_q
        break
      }
      q <- new_q
    }
    w_path[t, ] <- q
    
    if(adaptive_nu) {
      Q <- Q + (adaptive_steps[t] - 1) / adaptive_steps[t + 1] * (Q - old_Q)
    }
    
    if(progress) {
      setTxtProgressBar(pb, t)
    }
  }
  
  end <- Sys.time()
  
  if(progress) {
    close(pb)
  }
  
  return(
    list(
      U = Q,
      w = q,
      lambda = lambda,
      cobra_diffs = cobra_diffs,
      biconvex_diffs = biconvex_diffs,
      w_path = w_path,
      objective = objective_vals,
      rss = rss_vals,
      time = as.numeric(end - start)
    )
  )
}

mat_df <- function(bcbc_result,
                   cluster_w_weights = TRUE,
                   filter_weight = -Inf) {
  # filter weight will be as if columns with weight below filter weight did not exist
  w_prime_idx <- which(bcbc_result$w > filter_weight)
  w_prime <- bcbc_result$w[w_prime_idx]
  U_prime <- bcbc_result$U[, w_prime_idx]
  df <- as.data.frame(U_prime) |>
    mutate(row_num = row_number()) |>
    pivot_longer(!row_num, names_to = "col_num") |>
    mutate(col_num = as.numeric(str_replace(col_num, "V", "")),
           weight = rep(w_prime, max(row_num)))
  
  lambda <- bcbc_result$lambda
  # TODO should row dist be weighted?
  if (cluster_w_weights) {
    w_vals = sqrt(w_prime + lambda * w_prime)
    solution_mat = sweep(U_prime, 2, w_vals, "*")
  } else {
    solution_mat <- U_prime
  }
  row_dist <- dist(solution_mat)
  col_dist <- dist(t(solution_mat))
  
  h_row <- hclust(row_dist)
  h_col <- hclust(col_dist)
  df <- df |>
    inner_join(tibble(
      row_num = h_row$order,
      order_row = 1:length(h_row$order)
    )) |>
    inner_join(tibble(
      col_num = h_col$order,
      order_col = 1:length(h_col$order)
    ))
  
  df
}

centroid_rows <- function(mat, mat_for_dist, threshold, calculate_centroids=TRUE) {
  to_return <- matrix(NA, nrow=nrow(mat), ncol=ncol(mat))
  row_clusters <- dbscan::comps(frNN(mat_for_dist, threshold, sort=F))
  names(row_clusters) <- 1:nrow(mat)
  sorted_membership <- sort(row_clusters)
  node_indices <- as.integer(names(sorted_membership))
  cluster_sizes <- table(row_clusters)
  curr_index <- 1
  
  if(!calculate_centroids) {
    return(list(mat = to_return, cluster_info = row_clusters))
  }
  
  for (cluster_id in names(cluster_sizes)) {
    cluster_size <- cluster_sizes[cluster_id]
    cluster_members <-
      node_indices[curr_index:(curr_index + cluster_size - 1)]
    if (cluster_size > 1) {
      centroid <- colSums(mat[cluster_members,]) / cluster_size
      to_return[cluster_members,] <-
        rep(centroid, each = cluster_size)
    }
    
    curr_index <- curr_index + cluster_size
  }
  
  return(list(mat = to_return, cluster_info = row_clusters))
}

thresholded_solution <- function(bcbc_result,
                                 percent_of_noise,
                                 cluster_w_weights = TRUE) {
  to_return <- list()
  w <- bcbc_result$w
  U <- bcbc_result$U
  lambda <- bcbc_result$lambda
  if (cluster_w_weights) {
    w_vals <- sqrt(w + lambda * w)
    solution_mat <- sweep(U, 2, w_vals, "*")
  } else {
    solution_mat <- U
  }
  row_sd <- sd(sqrt(rowSums(solution_mat ^ 2)))
  threshold <- row_sd * percent_of_noise
  clustering <- centroid_rows(U, solution_mat, threshold)
  bcbc_result$U <- clustering$mat
  bcbc_result$cluster_info <- clustering$cluster_info
  
  to_return$w <- w
  to_return$lambda <- lambda
  to_return$cluster_info <- clustering$cluster_info
  to_return$U <- clustering$mat
  
  return(to_return)
}

get_row_clusters <- function(mat, threshold) {
  clustering <- centroid_rows(mat, mat, threshold, calculate_centroids = FALSE)
  
  return(clustering$cluster_info)
}

plot_matrix <- function(from_mat_df,
                        fill_attr = "value",
                        alpha_weight = FALSE,
                        bin_scale = FALSE) {
  if (alpha_weight) {
    raster <-
      geom_raster(aes(fill = .data[[fill_attr]], alpha = weight))
  } else {
    raster <- geom_raster(aes(fill = .data[[fill_attr]]))
  }
  
  if (!bin_scale) {
    plot_scale = scale_fill_gradient2()
  } else {
    plot_scale = scale_fill_steps2(n.breaks = bin_scale)
  }
  
  ggplot(from_mat_df, aes(x = order_col, y = order_row)) +
    raster +
    plot_scale +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

get_cv_metrics <- function(X,
                           bcbc_run,
                           weighted_clusters = TRUE,
                           percent_noise = 0.25) {
  U <- bcbc_run$U
  w <- bcbc_run$w
  lambda <- bcbc_run$lambda
  n <- nrow(X)
  p <- ncol(X)
  swept <- sweep(U, 2, w ^ 2 + lambda * w, "*")
  rss <- sum(sweep(X - U, 2, w ^ 2 + lambda * w, "*") ^ 2)
  rss_no_lambda_sq <- sum(sweep(X - U, 2, w ^ 2, "*") ^ 2)
  if (weighted_clusters) {  # TODO make sd calculate only using nonzero features
    row_sd <- sd(sqrt(rowSums(swept ^ 2)))
    col_sd <- sd(sqrt(colSums(swept ^ 2)))
    solution_mat <- swept
  } else {
    row_sd <- sd(sqrt(rowSums(U ^ 2)))
    col_sd <- sd(sqrt(colSums(U ^ 2)))
    solution_mat <- U
  }
  
  row_clusters <- get_row_clusters(solution_mat, row_sd * percent_noise)
  col_clusters <- get_row_clusters(t(solution_mat), col_sd * percent_noise)
  
  num_row_clusters <- length(unique(row_clusters))
  num_col_clusters <- length(unique(col_clusters))
  
  sparsity_kurtosis <- sum(w ^ 2) ^ 2 / sum(w ^ 4)
  sparsity_non_zero_w <- sum(w > 1e-10)
  sparsity_non_zero_sol <- sum(abs(swept) > 1e-10)
  return(list(
    cv_metrics = list(
      num_row_clusters = num_row_clusters,
      num_col_clusters = num_col_clusters,
      sparsity_kurtosis = sparsity_kurtosis,
      sparsity_non_zero_w = sparsity_non_zero_w,
      sparsity_non_zero_sol = sparsity_non_zero_sol,
      rss = rss,
      rss_no_lambda_sq = rss_no_lambda_sq,
      BIC = n * p * log(rss_no_lambda_sq / n * p) +
        1 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC1 = n * p * log(rss_no_lambda_sq / n * p) +
        2 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w),
      eBIC2 = n * p * log(rss_no_lambda_sq / n * p) +
        3 * log(n * p) * (num_row_clusters * num_col_clusters +
                            sparsity_non_zero_w)
    ),
    row_clusters = row_clusters,
    col_clusters = col_clusters
  ))
}

tune_bcbc <- function(X,
                      model = rscobra,
                      lambdas = c(1),
                      k_rows = c(2),
                      k_cols = c(2),
                      nus = c(1),
                      gammas = c(1),
                      tmaxs = c(NA),
                      tmax_outers = c(100),
                      tmax_cobras = c(100),
                      tmax_biconvexs = c(100),
                      phis = c(0.5),
                      recalculate_weights = c(TRUE),
                      tols = c(1e-6),
                      weighted_clusters = TRUE,
                      percent_noise = c(0.1),
                      progress = FALSE) {
  
  all_params <-
    expand.grid(
      lambda = lambdas,
      k_row = k_rows,
      k_col = k_cols,
      nu = nus,
      gamma = gammas,
      tmax = tmaxs,
      tmax_outer = tmax_outers,
      tmax_cobra = tmax_cobras,
      tmax_biconvex = tmax_biconvexs,
      phi = phis,
      recalculate_weights = recalculate_weights,
      tol = tols
    )
  
  cv_data <- data.frame(all_params) |>
    mutate(index=row_number())
  
  all_runs <- list()
  row_clusters <- list()
  col_clusters <- list()
  
  future_results <- future_lapply(
    1:nrow(all_params),
    future.seed = TRUE,
    function(param_set) {
      params <- all_params[param_set,]
      result <- do.call(model,
                        c(list(X=X, progress=progress), params))
      
      cv_metrics <- get_cv_metrics(X,
                                   result,
                                   weighted_clusters = weighted_clusters,
                                   percent_noise = percent_noise)
      
      list(
        result = result,
        cv_metrics = cv_metrics$cv_metrics,
        row_clusters = cv_metrics$row_clusters,
        col_clusters = cv_metrics$col_clusters
      )
    })
  
  for(param_set in 1:nrow(all_params)) {
    all_runs[[param_set]] = future_results[[param_set]]$result
    cv_metrics <- future_results[[param_set]]$cv_metrics
    cv_data[param_set, names(cv_metrics)] <- cv_metrics
    row_clusters[[param_set]] = future_results[[param_set]]$row_clusters
    col_clusters[[param_set]] = future_results[[param_set]]$col_clusters
  }
  
  return(list(all_runs = all_runs, cv_data = cv_data,
              row_clusters = row_clusters, col_clusters = col_clusters))
}

melt_cv_data <- function(tuning_data) {
  # Makes a big df mostly for visualization
  num_runs <- nrow(tuning_data$cv_data)
  to_return <-
    cbind(mat_df(tuning_data$all_runs[[1]]), tuning_data$cv_data[1, ])
  for (i in 2:num_runs) {
    to_return <- rbind(to_return,
                       cbind(mat_df(tuning_data$all_runs[[i]]), tuning_data$cv_data[i, ]))
  }
  return(to_return)
}

gen_checkerboard <- function(n,
                             p,
                             num_row_clusters,
                             num_col_clusters,
                             p_extra = 0,
                             noise = 1,
                             shuffle = TRUE) {
  row_partition <- sort(sample(1:num_row_clusters, n, replace = TRUE))
  col_partition <-
    sort(sample(1:num_col_clusters, p, replace = TRUE))
  mu_kr <- matrix(
    runif(num_row_clusters * num_col_clusters,-2, 2),
    nrow = num_row_clusters,
    ncol = num_col_clusters
  )
  
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
  if (p_extra > 0) {
    data_mat[, (p + 1):(p + p_extra)] <-
      rnorm(n * p_extra, mean = 0, sd = noise)
    centers[, (p + 1):(p + p_extra)] <- 0
    col_partition <- c(col_partition, rep(num_col_clusters + 1, p_extra))
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
    return(list(X = scale(shuffled_row), centers = shuffled_center, 
                row_partition=row_partition, col_partition=col_partition))
  }
  return(list(X = scale(data_mat), centers = centers,
              row_partition=row_partition, col_partition=col_partition))
}

# From https://github.com/kharchenkolab/vrnmf/tree/main
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

palm <- function(X,
                 lambda,
                 k_row = 4,
                 k_col = 4,
                 nu = 1,
                 gamma = 10,
                 phi = 0.05,
                 tmax = NA, # Set this to set both below
                 tmax_cobra = 100,
                 tmax_biconvex = NA, # Unused, here for compatibility with rscobra
                 tmax_outer = 100,
                 tol = 1e-6,
                 recalculate_weights = TRUE,
                 approx = 0,
                 threshold = 0,
                 progress = TRUE) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  U <- X
  w <- rep(1/p, p)
  
  if (!is.na(tmax)) {
    tmax_cobra <- tmax
    tmax_outer <- tmax
  }
  
  cobra_diffs <- rep(NA, tmax_outer)
  w_diffs <- rep(NA, nrow = tmax_outer)
  w_path <- matrix(NA, nrow = tmax_outer, ncol = p)
  objective_vals <- rep(NA, tmax_outer)
  rss_vals <- rep(NA, tmax_outer)
  
  if (progress) {
    pb <- txtProgressBar(min = 0,
                         max = tmax_outer,
                         initial = 0)
  }
  
  start <- Sys.time()
  for (t in 1:tmax_outer) {
    nu_U = nu
    U_step = U - nu_U * sweep(X - U, 2, w ^ 2 + lambda * w, "*")
    if (recalculate_weights || t == 1) {
      wts <- fast_gkn_weights(
        t(U_step),
        k_row = k_row,
        k_col = k_col,
        phi = phi,
        approx = approx
      )
      w_row <- wts$w_row
      w_col <- wts$w_col
      E_row <- wts$E_row
      E_col <- wts$E_col
      row_fusion <- wts$row_fusion
      col_fusion <- wts$col_fusion
    }
    
    rss_vals[t] <- sum(sweep(X - U, 2, w ^ 2 + lambda * w, "*") ^ 2)
    objective_vals[t] <-
      gamma * (row_fusion + col_fusion) + rss_vals[t] / 2
    
    # Cobra has rows as features, columns as samples
    cobra_result <-
      cobra(t(U_step),
            E_row,
            E_col,
            w_row,
            w_col,
            gamma = gamma,
            max_iter = tmax_cobra,
            tol=tol)
    
    U_prime <- t(cobra_result$U[[1]])
    U_diff <- sum(abs(U_prime - U)) / sum(abs(U_prime))
    cobra_diffs[t] <- U_diff
    if (U_diff < tol) {
      U <- U_prime
      break
    }
    U <- U_prime
    
    col_sum_sq <- colSums((X - U) ^ 2)
    lipschitz_fixed_U = sqrt(sum(col_sum_sq ^ 2))
    nu_w = 1/(1.1 * lipschitz_fixed_U)
    
    w_prime <-
      projection_onto_simplex(w - nu_w * (w + lambda / 2) * col_sum_sq)
    
    w_diff <- sum((w_prime - w) ^ 2)
    w_diffs[t] <- w_diff
    w <- w_prime
    w_path[t, ] <- w
    if (progress) {
      setTxtProgressBar(pb, t)
    }
  }
  
  end <- Sys.time()
  
  if (progress) {
    close(pb)
  }
  
  return(
    list(
      U = U,
      w = w,
      lambda = lambda,
      cobra_diffs = cobra_diffs,
      w_diffs = w_diffs,
      w_path = w_path,
      objective = objective_vals,
      rss = rss_vals,
      time = as.numeric(end - start)
    )
  )
}