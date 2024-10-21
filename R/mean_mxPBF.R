mxPBF_mean <- function(given_data, nws = c(25, 60, 100), alps = seq(0, 10, 0.05), FPR_want = 0.05, n_sample = 300, n_cores = 1) {
  results <- lapply(nws, function(nw) {
    cp <- numeric(0)
    mxPBF_simulated <- simulate_mxPBF_mean2(given_data, nw, alps, n_sample, n_cores)
    alp_selected <- alps[which.min(abs(rowSums(exp(mxPBF_simulated) > 10) / n_sample - FPR_want))]
    bf_given_data <- cpd_mean_mxPBF2(given_data, nw, alp_selected, n_cores)
    exp_bf <- exp(bf_given_data)
    while (any(exp_bf > 10)) {
      i_tilde <- which(exp_bf > 10)[1]
      i_hat <- i_tilde + which.max(bf_given_data[i_tilde:(i_tilde + nw - 1)]) - 1
      cp <- c(cp, i_hat + nw)
      exp_bf[i_tilde:(i_hat + nw - 1)] <- 0
    }
    mxPBF_result <- list("Change_points" = cp,
                         "Bayes_Factors" = bf_given_data,
                         "Selected_alpha" = alp_selected,
                         "Window_size" = nw)
    class(mxPBF_result) <- "mxPBF"
    return(mxPBF_result)
  })
  names(results) <- paste("Window_size", nws, sep = "_")
  return(results)
}
