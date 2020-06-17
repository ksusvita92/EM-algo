# EM-algo-

The folder project consist of:

R scripts:

  + **PSA_EM algo.R**: stands for Pairwise Survival Analysis EM-algorithm. An R script to estimate gamma parameter (i.e. shape and scale) using EM-algorithm. We use MoM (Method of Moment) to estimate the parameters. Required inputs are a vector of data points (representing serial interval) and number of mixture components. The default setting is with tol=1e-5 and 1000 iterations. 
  + **PSA_clustering.R**: stands for Pairwise Survival Analysis clustering. An R script for clustering. The required inputs are the distance matrix and the threshold. Cases will be clustered if their distances are within the given threshold. The function returns a list of cluster members.
  + **PSA_time phylo.R**: A script to produce a collection of subtrees. Each subtree represents phylo tree of cases in the same cluster. If we do not provide a tree, we can use a metadata that represents a collection of serial intervals.
