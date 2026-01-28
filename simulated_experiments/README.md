Here we describe some of the implementations for all methods compared to BCBC. See `simulated_experiments.qmd` for more details.

-   Adaptive BCBC: Optimization of BCBC with adaptive weights using a normalized GKNN graph with $k=5$ and $\tau = 1$, updated every iteration for both rows and columns.
-   Approximate Adaptive Nearest Neighbor BCBC: Optimization of BCBC with adaptive weights using an approximate (via `RcppHNSW` on CRAN) normalized GKNN graph with $k=5$ and $\tau = 1$, updated every iteration for both rows and columns.
-   BCBC: Optimization of BCBC with PALM using a normalized GKNN graph with $k=5$ and $\tau = 1$ for both the rows and columns.
-   BCEL (Available at [yanzhong07/BCEL](https://github.com/yanzhong07/BCEL)): Optimization of a biconvex objective with exclusive LASSO regularization. As implemented by the authors, selection of the regularization hyperparameters is done via stability selection.
-   COBRA (Archived on CRAN as `cvxbiclustr`): Optimization of convex biclustering using a normalized GKNN graph with $k=5, \tau = 1$. Tuning of $\gamma$ is performed with a hold-out set as recommended in \citet{chiConvexBiclustering2017}, with final cluster memberships determined with a thresholding radius of 0.25 times the standard deviation of pairwise row and column distances. To match the use of COBRA in practice, we optimize on only the top 200 (the true number of informative) features by variance.
-   SC-Biclust (Archived on CRAN as `SCBiclust`): Repeatedly performs the sparse clustering method of "A framework for feature selection in clustering" (Witten and Tibshirani, 2010) and uses a hypothesis test to filter features to a bicluster. Stability selection is used to determine the number of biclusters.
-   SparseBC (Archived on CRAN as `sparseBC`): A generalization of $k$-means for biclustering with a LASSO penalty to promote sparsity. The number of biclusters is calculated via a hold-out set as described in the original paper. The regularization parameter for the LASSO penalty was chosen via BIC.

In addition, simulations exist for the following methods, which were not included due to poor performance:
-   Dynamic Tree Cut (Archived on CRAN as `dynamicTreeCut`): Complete-linkage hierarchical clustering which is adaptively pruned into clusters. This method is applied independently to the rows and columns with the true minimum cluster size provided as a hyperparameter.
-   PLAID (Archived on CRAN as `biclust`): Heuristic search method combined with additive model as implemented in the `biclust` package.
-   S4VD (Archived on CRAN as `s4vd`): Sparse singular value decomposition with stability selection. All hyperparameters are provided to be as informative as possible such as cluster size and number of biclusters.
