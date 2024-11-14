Dette_cp <- function(data) {
  library(Matrix)

  Data2Tmat <- function(data) {
    n <- nrow(data)
    p <- ncol(data)
    Tmat <- matrix(0, p * (p + 1) / 2, n)
    for (k in 1:n) {
      row <- data[k, ]
      outer <- row %*% t(row)
      lower_indices <- lower.tri(outer, diag = TRUE)
      Tmat[, k] <- outer[lower_indices]
    }
    return(Tmat)
  }

  Tmat2Dvec <- function(Tmat) {
    p2 <- nrow(Tmat)
    n <- ncol(Tmat)
    Tsum <- rowSums(Tmat)
    left_cumsum1 <- t(apply(Tmat, 1, cumsum))
    left_cumsum2 <- left_cumsum1^2
    left_cumsum3 <- t(apply(Tmat^2, 1, cumsum))

    right_cumsum1 <- matrix(rep(Tsum, each=n), nrow = p2, byrow = TRUE) - left_cumsum1
    right_cumsum2 <- right_cumsum1^2
    right_cumsum3 <- matrix(rep(rowSums(Tmat^2), each=n), nrow = p2, byrow = TRUE) - left_cumsum3

    mixed_cumsum <- left_cumsum1 * right_cumsum1

    a <- numeric(n)
    b <- numeric(n)
    c <- numeric(n)

    for (k in 2:(n-2)) {
      a[k] = (n - k) / (n * (k - 1) * (n - 3))
      b[k] = k / (n * (n - k - 1) * (n - 3))
      c[k] = 2 / (n * (n - 3))
    }

    a1 <- matrix(rep(a, each = p2), nrow = p2)
    b1 <- matrix(rep(b, each = p2), nrow = p2)
    c1 <- matrix(rep(c, each = p2), nrow = p2)

    A <- a1 * (left_cumsum2 - left_cumsum3)
    A <- A[, -c(1, n-1, n)]
    T1 <- rowSums(A)

    B <- b1 * (right_cumsum2 - right_cumsum3)
    B <- B[, -c(1, n-1, n)]
    T2 <- rowSums(B)

    C <- c1 * mixed_cumsum
    C <- C[, -c(1, n-1, n)]
    T3 <- rowSums(C)

    result <- T1 + T2 - T3

    return(result)
  }

  select_threshold <- function(data) {
    n <- nrow(data)
    p <- ncol(data)
    p2 <- p * (p + 1) / 2

    data_dot <- data - colMeans(data)
    Tmat <- Data2Tmat(data_dot)
    rm(data_dot) # Remove unused variable
    D <- Tmat2Dvec(Tmat)

    Z <- matrix(0, ncol = p2, nrow = floor(n/2))
    for (k in 1:floor(n/2)) {
      Z[k ,] <- Tmat[, 2 * k] - Tmat[, 2 * k - 1]
    }
    rm(Tmat) # Remove unused variable

    O <- 1 / sqrt(2) * apply(Z, 2, sd)
    rm(Z) # Remove unused variable

    Y <- matrix(rnorm(p2 * n), nrow = p2, ncol = n)
    X_star <- sweep(Y, 1, O, `*`)
    rm(Y, O) # Remove unused variables

    D_star <- Tmat2Dvec(X_star)
    rm(X_star) # Remove unused variable
    gc() # Garbage collection for large variable

    tau <- max(D_star)
    loc <- which(D > tau)
    len <- length(loc)

    return(list(location = loc, length = len))
  }

  p <- ncol(data)
  n <- nrow(data)

  threshold <- select_threshold(data)
  loc <- threshold$location
  len <- threshold$length
  if (len == 0) {
    return(integer(0))
  }
  data_dot <- data - colMeans(data)
  Tmat <- Data2Tmat(data_dot)
  rm(data_dot) # Remove unused variable
  T_reduced <- Tmat[loc, , drop=FALSE]
  rm(Tmat) # Remove unused variable
  gc() # Garbage collection for large variable

  U <- numeric(n)
  for (k in 2:(n-2)) {
    T_left <- T_reduced[, 1:k, drop=FALSE]
    T_right <- T_reduced[, (k+1):n, drop=FALSE]

    s1 <- rowSums(T_left)
    Q1 <- (1 / (k * (k - 1))) * (sum(s1^2) - sum(T_left^2))

    s2 <- rowSums(T_right)
    Q2 <- (1 / ((n - k) * (n - k - 1))) * (sum(s2^2) - sum(T_right^2))

    Q3 <- (1 / (k * (n - k))) * sum(s1 * s2)
    Q <- Q1 + Q2 - 2 * Q3

    U[k] <- (k * (k - 1) * (n - k) * (n - k - 1)) / n^4 * Q
  }
  k_hat <- which.max(U)

  return(k_hat)
}
