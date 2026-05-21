# BCBC

Code for Biconvex Biclustering from https://arxiv.org/abs/2604.03936

## Installation

```{r}
devtools::install_github("SamGRosen/BCBC") 
```

## Complete Example

We show how to generate some data, do proper tuning according to the paper, and examine the results.


### Simulate Data

We simulate some data with a checkerboard structure, but with many uninformative features.

```{r}
n = 100
p = 150 + 300

checker <- gen_checkerboard(n, p - 300, 5, 5, p_extra = 300, shuffle = TRUE)
as_df <- matrix_fit_to_df(checker$X, labels = checker$row_partition)
plot_fit(as_df)

```

### Gamma Tuning

Following the tuning in the paper, we first tune gamma with a hold-out set.

```{r}
gamma_cv <- cv.BCBC_holdout(
  checker$X,
  holdout_size = 0.25,
  lambdas = 0,
  gammas = n * 1.5 ^ seq(-4, 4, 1),
  k_samples = 4,
  k_features = 4,
  phis = 1,
  tols = 10 ^ -3,
  tmax_hierarchy = c(3, 50, 100),
  recalculate_weights = c(TRUE)
)

best_gamma <- gamma_cv$gamma_cv_data |> slice(which.min(rss_heldout))
```

### Lambda Tuning

Afterwards, with the optimal gamma, we now perform tuning of lambda using eBIC.

```{r}
lambda_cv <- cv.BCBC(
  checker$X,
  gamma = c(best_gamma$gamma),
  lambdas = seq(0, 120, 8) / p,
  k_samples = 4,
  k_features = 4,
  phi = 1,
  percent_noise = seq(0.025, 0.25, 0.025),
  recalculate_weights = TRUE,
  tols = 1e-4,
  tmax_outer = 50,
  tmax_cobra = 100
)

```

### Displaying Results

Several utilities are available to plot results.

```{r}

best_run_info <- lambda_cv$cv_data |>
   filter(is.finite(eBIC2)) |>
   slice(which.min(eBIC2))

best_param_index <- best_run_info$param_index
best_cluster_index <- best_run_info$index

row_groups <- lambda_cv$row_clusters[[best_cluster_index]]
col_groups <- lambda_cv$col_clusters[[best_cluster_index]]

best_run <- lambda_cv$all_runs[[best_param_index]]

fit_as_df <- bcbc_result_to_df(best_run, filter_weight = 0)
plot_fit(fit_as_df)
```


## Run Simulations

Simulations can be run with `BCBC/simulated_experiments/simulated_experiments.qmd`
or `BCBC/simulated_experiments/run_all_simulations.sh`.

## Run Real Data Example

Real data can be run with `BCBC/real_data_experiments/run_all_experiments.sh`.

The lymphoma experiment can be run with `BCBC/real_data_experiments/lymphoma_experiment.sh`.

## Other methods

To run the simulated experiments, the following packages need to be installed.

```{r}
install.packages("SCBiclust")
install.packages("sparseBC")
install.packages("biclust")
install.packages("s4vd")
install.packages("dynamicTreeCut")

devtools::install_github("yanzhong07/BCEL") 
devtools::install_github("cran/cvxbiclustr") 
```

## Citation

Please cite the following paper when using the code from this repository:

```
@misc{rosen2026biconvexbiclustering,
      title={Biconvex Biclustering}, 
      author={Sam Rosen and Eric C. Chi and Jason Xu},
      year={2026},
      eprint={2604.03936},
      archivePrefix={arXiv},
      primaryClass={stat.ML},
      url={https://arxiv.org/abs/2604.03936}, 
}
```
