rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
library(changepoint.geo)
library(changepoint.np)
library(InspectChangepoint)
library(ecp)
library(dplyr)
library(RSpectra)
library(hdbcp)
source("util_funcs.R")

# Parameters
## Simulation
n_rep <- 50 # The number of simulation replications
n <- 500 # The number of observations
pre_value <- 0.3 # The value of nonzero precision matrix
single_cp <- 250 # Location of the change point for single change point scenarios
multiple_cp <- c(150,300,350) # Locations of change points for multiple change point scenarios

## mxPBF
n_cores <- 1 # The number of cores for parallelization
n_sample <- 300 # The number of sample generations for selecting the hyperparameter alpha
fpr_want = 0.05 # The prespecified FPR for selecting the hyperparameter alpha
alps <- seq(1,10,0.05) # The grid of alpha values for selecting the hyperparameter alpha
nws <- c(25, 60, 100) # The set of window sizes

# Settings
signals <- c(0.1, 0.3, 0.5, 0.7, 1, 1.2, 1.5, 2) # The size of signals at change points
pre_proportions <- c(0.01,0.4) # The proportion of non-zero elements in the precision matrix
dims <- c(200, 500, 800) # The number of variables (Dimension)
settings <- expand.grid(signals, pre_proportions, dims)

# Data frame to store the results
methods <- c("mxPBF_major", "GeomCP", "Inspect", "Edivisive")
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
  pre_proportion <- settings[set,2]
  p <- settings[set,3]

  for (r in 1:n_rep) {
    given_datasets <- generate_mean_datasets(n, p, signal_size, pre_proportion, pre_value,
                                             single_point = single_cp, multiple_points = multiple_cp)
    for (t in 1:length(data_types_df)) {
      given_data <- matrix(given_datasets[,,t], n, p)

      ## mxPBF method
      res_mxPBF <- mxPBF_mean(given_data, nws, alps, fpr_want, n_sample, n_cores)

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
      cps_geom <- geomcp(given_data, msl = nws[1], penalty = "MBIC")
      cp_distance <- cpt.np(cps_geom@distance, penalty = "MBIC",
                            method="PELT", test.stat="empirical_distribution",
                            class=TRUE, minseglen=nws[1],
                            nquantiles =4*log(n))
      cp_angle <- cpt.np(cps_geom@angle, penalty = "MBIC",
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

      ## Inspect method
      out <- capture.output(cps_inspect <- inspect(t(given_data), threshold = compute.threshold(n, p, show_progress = F), show_progress = F))
      index_cps <- which(results_cps$Setting == set &
                           results_cps$Replication == r &
                           results_cps$Data_Type == data_types_df[t] &
                           results_cps$Method == "Inspect")
      index_others <- which(results_others$Setting == set &
                              results_others$Replication == r &
                              results_others$Data_Type == data_types_df[t] &
                              results_others$Method == "Inspect")
      results_cps$Detected_Change_Points[[index_cps]] <- inspect_prune(cps_inspect$changepoints, nws[1])
      results_others$info[[index_others]] <- I(list(cps_inspect$changepoints))


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
    }
    print(paste("Replication", r, "in Setting", set, "completed"))
  }
  save(results_cps, results_others, file = paste("result_mean_setting_to_",set,".RData",sep=""))
}
