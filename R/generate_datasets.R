# This function generates a simulated dataset containing change points in the mean vector.

generate_mean_datasets <- function(n = 500, p = 200, signal_size = 1, pre_proportion = 0.4, pre_value = 0.3, single_point = round(n/2), multiple_points = c(round(n/4), round(2*n/4), round(3*n/4)), type = c(1,2,3,4,5)) {
  generate_mean_data <- function(n, mu, mu2, cov, signal_loc) {
    signal_loc <- sort(signal_loc)
    start <- 1
    data_list <- list()
    for (i in seq_along(signal_loc)) {
      end <- signal_loc[i]
      if (i %% 2 == 1) {
        data_list[[i]] <- cpp_mvrnorm(end - start + 1, mu, sigma)
      } else {
        data_list[[i]] <- cpp_mvrnorm(end - start + 1, mu2, sigma)
      }
      start <- end + 1
    }
    if (start <= n) {
      if (length(signal_loc) %% 2 == 1) {
        data_list[[length(signal_loc) + 1]] <- cpp_mvrnorm(n - start + 1, mu2, sigma)
      } else {
        data_list[[length(signal_loc) + 1]] <- cpp_mvrnorm(n - start + 1, mu, sigma)
      }
    }
    data <- do.call(rbind, data_list)

    return(data)
  }
  # Initialize the mean vector
  mu <- rep(0, p)
  mu_R <- mu_M <- mu

  # Set the indices where the signal will be applied
  indices_R <- sample(length(mu), 5)
  indices_M <- sample(length(mu), p / 2)

  # Add signal to the mean vectors
  mu_R[indices_R] <- signal_size
  mu_M[indices_M] <- signal_size

  # Create the Omega matrix
  omega <- matrix(0, p, p)
  lower_tri_indices <- which(lower.tri(omega))
  selected_indices <- sample(lower_tri_indices, floor(pre_proportion * length(lower_tri_indices)))
  omega[selected_indices] <- pre_value
  omega <- omega + t(omega)
  diag(omega) <- diag(omega) / 2

  # Adjust the Omega matrix if the minimum eigenvalue is too small
  min_eigen_val <- min(eigen(omega)$values)
  if (min_eigen_val < 1e-5) {
    omega <- omega + (-min_eigen_val + 0.1^3) * diag(1, p)
  }

  # Calculate the Sigma matrix
  sigma <- solve(omega)

  # Generate data
  sim_data_H0 <- cpp_mvrnorm(n, mu, sigma)
  sim_data_H1R_single <- generate_mean_data(n, mu, mu_R, sigma, single_point)
  sim_data_H1M_single <- generate_mean_data(n, mu, mu_M, sigma, single_point)
  sim_data_H1R_multiple <- generate_mean_data(n, mu, mu_R, sigma, multiple_points)
  sim_data_H1M_multiple <- generate_mean_data(n, mu, mu_M, sigma, multiple_points)

  # Return the results as an array
  given_datasets <- array(c(sim_data_H0, sim_data_H1R_single, sim_data_H1M_single,
                            sim_data_H1R_multiple, sim_data_H1M_multiple),
                          dim = c(dim(sim_data_H0), 5))

  return(given_datasets[, , type])
}

# This function generates a simulated dataset containing change points in the covariance matrix.

generate_cov_datasets <- function(n = 500, p = 200, signal_size = 1, sparse = TRUE, single_point = round(n/2), multiple_points = c(round(n/4), round(2*n/4), round(3*n/4)), type = c(1,2,3,4,5)) {
  gen_U <- function(p, signal, rare = TRUE) {
    if (rare) {
      U <- matrix(0, nrow = p, ncol = p)

      lower_tri_indices <- which(lower.tri(U ,diag = T))
      selected_indices <- sample(lower_tri_indices, 5)

      U[selected_indices] <- runif(5, 0, signal)
      U <- U + t(U)
      diag(U) <- diag(U) / 2
    } else {
      u <- runif(p, 0, signal)
      U <- u %*% t(u)
    }
    return(U)
  }
  gen_Sigma <- function(p, sparse = TRUE) {
    if (sparse) {
      Delta_1 <- matrix(0, nrow = p, ncol = p)
      lower_tri_indices <- which(lower.tri(Delta_1 ,diag = T))
      selected_indices <- sample(lower_tri_indices, floor(0.05 * length(lower_tri_indices)))
      Delta_1[selected_indices] <- 0.5
      Delta_1 <- Delta_1 + t(Delta_1)
      diag(Delta_1) <- diag(Delta_1) / 2
      min_eigenvalue <- min(eigen(Delta_1)$values)
      if (min_eigenvalue <= 1e-5) {
        Delta_1 <- Delta_1 + (abs(min_eigenvalue) + 0.05) * diag(p)
      }
      d <- runif(p, 0.5, 2.5)
      D_sqrt <- diag(sqrt(d))
      Sigma_1 <- D_sqrt %*% Delta_1 %*% D_sqrt
    } else {
      omega <- runif(p, 1, 5)
      Omega <- diag(omega)
      Delta <- matrix(0, nrow = p, ncol = p)
      for (i in 1:p) {
        for (j in 1:p) {
          Delta[i, j] <- (-1)^(i + j) * 0.4^(abs(i - j)^(1/10))
        }
      }
      Sigma_1 <- Omega %*% Delta %*% Omega
    }
    return(Sigma_1)
  }
  generate_cov_data <- function(n, mu, sigma, sigma2, signal_loc) {
    signal_loc <- sort(signal_loc)
    start <- 1
    data_list <- list()
    for (i in seq_along(signal_loc)) {
      end <- signal_loc[i]
      if (i %% 2 == 1) {
        data_list[[i]] <- cpp_mvrnorm(end - start + 1, mu, sigma)
      } else {
        data_list[[i]] <- cpp_mvrnorm(end - start + 1, mu, sigma2)
      }
      start <- end + 1
    }
    if (start <= n) {
      if (length(signal_loc) %% 2 == 1) {
        data_list[[length(signal_loc) + 1]] <- cpp_mvrnorm(n - start + 1, mu, sigma2)
      } else {
        data_list[[length(signal_loc) + 1]] <- cpp_mvrnorm(n - start + 1, mu, sigma)
      }
    }
    data <- do.call(rbind, data_list)

    return(data)
  }
  mu <- rep(0,p)
  sigma <- gen_Sigma(p,sparse)
  U_rare <- gen_U(p,signal_size,rare=T)
  U_many <- gen_U(p,signal_size,rare=F)
  sigma_R <- sigma + U_rare
  sigma_M <- sigma + U_many
  min_eigen1 <- sigma %>% eigen %>% .$values %>% min
  min_eigen2 <- sigma_R %>% eigen %>% .$values %>% min
  min_eigen3 <- sigma_M %>% eigen %>% .$values %>% min
  if (any(c(min_eigen1,min_eigen2,min_eigen3) < 1e-5)) {
    delta <- abs(min(c(min_eigen1,min_eigen2,min_eigen3))) + 0.05
    sigma <- sigma + delta * diag(p)
    sigma_R <- sigma_R + delta * diag(p)
    sigma_M <- sigma_M + delta * diag(p)
  }

  # Generate data
  sim_data_H0 <- cpp_mvrnorm(n, mu, sigma)
  sim_data_H1R_single <- generate_cov_data(n, mu, sigma, sigma_R, single_point)
  sim_data_H1M_single <- generate_cov_data(n, mu, sigma, sigma_M, single_point)
  sim_data_H1R_multiple <- generate_cov_data(n, mu, sigma, sigma_R, multiple_points)
  sim_data_H1M_multiple <- generate_cov_data(n, mu, sigma, sigma_M, multiple_points)

  # Return the results as an array
  given_datasets <- array(c(sim_data_H0, sim_data_H1R_single, sim_data_H1M_single,
                            sim_data_H1R_multiple, sim_data_H1M_multiple),
                          dim = c(dim(sim_data_H0), 5))

  return(given_datasets[, , type])
}
