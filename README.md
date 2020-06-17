# EM-algo-

The folder project consist of:

R scripts:

  + **PSA_EM algo.R**: stands for Pairwise Survival Analysis EM-algorithm. An R script to estimate gamma parameter (i.e. shape and scale) using EM-algorithm. We use MoM (Method of Moment) to estimate the parameters. Required inputs are a vector of data points (representing serial interval) and number of mixture components. The default setting is with tol=1e-5 and 1000 iterations. 

