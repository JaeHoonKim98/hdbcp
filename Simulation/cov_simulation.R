rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(changepoint.geo)
library(changepoint.np)
library(ecp)
library(dplyr)
library(hdbcp)
source("cov_dette.R")
source("util_funcs.R")

# Parameters
## Simulation
n_rep <- 50 # The number of simulation replications
n <- 500 # The number of observations
single_cp <- 250 # Location of the change point for single change point scenarios
multiple_cp <- c(150,300,350) # Locations of change points for multiple change point scenarios

## mxPBF
n_cores <- 1 # The number of cores for parallelization
n_sample <- 300 # The number of sample generations for selecting the hyperparameter alpha
fpr_want = 0.05 # The prespecified FPR for selecting the hyperparameter alpha
alps <- seq(1,10,0.05) # The grid of alpha values for selecting the hyperparameter alpha
a0 <- 0.01 # The hyperparameter
b0 <- 0.01 # The hyperparameter
nws <- c(25, 60, 100) # The set of window sizes

# Settings
signals <- c(0.5, 1, 1.3, 2, 4, 6, 8, 10) # The size of signals at change points
sparses <- c(TRUE, FALSE) # Sparsity level of the covariance matrix
dims <- c(200, 500, 800) # The number of variables (Dimension)
settings <- expand.grid(signals, sparses, dims)

# Data frame to store the results
methods <- c("mxPBF_major", "GeomCP", "Dette", "Edivisive")
data_types_df <- c("H0", "H1_rare_single", "H1_many_single", "H1_rare_multiple", "H1_many_multiple")
settings_df <- 1:nrow(settings)
replications_df <- 1:n_rep

## Data frame containing detected change points
results_cps <- expand.grid(
  Setting = settings_df,
  Replication = replications_df,
  Data_Type = data_types_df,
  Method = methods,
  stringsAsFactors = FALSE
)

## Data frame containing infomations
results_others <- expand.grid(
  Setting = settings_df,
  Replication = replications_df,
  Data_Type = data_types_df,
  Method = methods,
  stringsAsFactors = FALSE
)
results_cps$Detected_Change_Points <- vector("list", length = nrow(results_cps))
results_others$info <- vector("list", nrow(results_others))


# Simulation
for (set in 1:nrow(settings)) {
  signal_size <- settings[set,1]
  sparse <- settings[set,2]
  p <- settings[set,3]

  for (r in 1:n_rep) {
    given_datasets <- generate_cov_datasets(n, p, signal_size, sparse, single_point = 250, multiple_points = c(150,300,350))

    for (t in 1:length(data_types_df)) {
      given_data <- matrix(given_datasets[,,t], n, p)

      ## mxPBF method
      res_mxPBF <- mxPBF_cov(given_data = given_data, a0 = a0, b0 = b0, nws = nws, alps = alps, FPR_want = fpr_want, n_sample = n_sample, n_cores = n_cores, centering = "skip")

      index_cps <- which(results_cps$Setting == set &
                           results_cps$Replication == r &
                           results_cps$Data_Type == data_types_df[t] &
                           results_cps$Method == "mxPBF_major")
      index_others <- which(results_others$Setting == set &
                              results_others$Replication == r &
                              results_others$Data_Type == data_types_df[t] &
                              results_others$Method == "mxPBF_major")
      results_cps$Detected_Change_Points[[index_cps]] <- majority_rule_mxPBF(res_mxPBF)
      results_others$info[[index_others]] <- I(res_mxPBF)


      ## GeomPELT method
      cp_distance <- cpt.np(distance.mapping(given_data), penalty = "MBIC",
                            method="PELT", test.stat="empirical_distribution",
                            class=TRUE, minseglen=nws[1],
                            nquantiles =4*log(n))
      cp_angle <- cpt.np(angle.mapping(given_data), penalty = "MBIC",
                         method="PELT", test.stat="empirical_distribution",
                         class=TRUE, minseglen=nws[1],
                         nquantiles =4*log(n))
      index_cps <- which(results_cps$Setting == set &
                           results_cps$Replication == r &
                           results_cps$Data_Type == data_types_df[t] &
                           results_cps$Method == "GeomCP")
      index_others <- which(results_others$Setting == set &
                              results_others$Replication == r &
                              results_others$Data_Type == data_types_df[t] &
                              results_others$Method == "GeomCP")
      results_cps$Detected_Change_Points[[index_cps]] <- geomcp_reconcile(cp_distance@cpts[-length(cp_distance@cpts)], cp_angle@cpts[-length(cp_angle@cpts)])
      results_others$info[[index_others]] <- I(list("distance" = cps_geom@distance, "angle" = cps_geom@angle))


      ## Edivisive method
      cps_edv <- e.divisive(given_data, R = 499, min.size = nws[1], alpha = 1)
      index_cps <- which(results_cps$Setting == set &
                           results_cps$Replication == r &
                           results_cps$Data_Type == data_types_df[t] &
                           results_cps$Method == "Edivisive")
      index_others <- which(results_others$Setting == set &
                              results_others$Replication == r &
                              results_others$Data_Type == data_types_df[t] &
                              results_others$Method == "Edivisive")
      if (length(cps_edv$estimates) > 2) {
        results_cps$Detected_Change_Points[[index_cps]] <- cps_edv$estimates[2:(length(cps_edv$estimates) - 1)]
      } else {
        results_cps$Detected_Change_Points[[index_cps]] <- integer(0)
      }
      results_others$info[[index_others]] <- I(cps_edv)

      ## Dette method (skip for mutiple change points scenarios)
      if (t != 4 && t != 5) {
        cps_Dette <- Dette_cp(given_data)
        index_cps <- which(results_cps$Setting == set &
                             results_cps$Replication == r &
                             results_cps$Data_Type == data_types_df[t] &
                             results_cps$Method == "Dette")
        results_cps$Detected_Change_Points[[index_cps]] <- cps_Dette
      }
    }
    print(paste("Replication", r, "in Setting", set, "completed"))
  }
  save(results_cps, results_others, file = paste("result_cov_setting_to_",set,".RData",sep=""))
}
