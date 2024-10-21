mxPBF_cov <- function(given_data, a0 = 0.01, b0 = 0.01, nws = c(25, 60, 100), alps = seq(0, 10, 0.05), FPR_want = 0.05, n_sample = 300, n_cores = 1, centering = "skip"){
  if (centering == "mean") {
    means <- colMeans(given_data)
    centered_data <- sweep(given_data, 2, means, FUN = "-")
  }
  if (centering == "median") {
    library(matrixStats)
    medians <- colMedians(given_data)
    centered_data <- sweep(given_data, 2, medians, FUN = "-")
  }
  if (centering == "skip") {
    centered_data <- given_data
  }

  results <- lapply(nws, function(nw) {
    cp <- numeric(0)
    mxPBF_simulated <- simulate_mxPBF_cov2(centered_data, a0, b0, nw, alps, n_sample, n_cores)
    alp_selected <- alps[which.min(abs(rowSums(exp(mxPBF_simulated) > 10) / n_sample - FPR_want))]
    bf_given_data <- cpd_cov_mxPBF2(centered_data, a0, b0, nw, alp_selected, n_cores)
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
