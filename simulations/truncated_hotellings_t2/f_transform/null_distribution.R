################################### M1

################################### Independent Error

sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45
p_val_star <- matrix(nrow = 10000, ncol = 6)

for (z in 1:6) {
  for (n in 1:10000) {
    t <- seq(from = 0, to = 1, length.out = 25)
    x <- rbind(
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z])
    )

    y <- rbind(
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma[z])
    )

    nx <- nrow(x)
    ny <- nrow(y)
    delta <- colMeans(x) - colMeans(y)
    p <- ncol(x)
    Sx <- cov(x)
    Sy <- cov(y)
    S_pooled <- ((nx - 1) * Sx + (ny - 1) * Sy) / (nx + ny - 2)

    L <- eigen(S_pooled)$values
    V <- eigen(S_pooled)$vectors
    A <- matrix(rep(0), nrow = 25, ncol = 25)

    val <- c()
    for (i in 1:25) {
      val[i] <- L[i] / L[1]
    }

    k <- tail(which(val > 10^-16), 1)

    t_val <- c()
    for (i in 1:k) {
      L_prime <- as.numeric(solve(L[i]))
      A1 <- L_prime * V[, i] %*% t(V[, i])
      A <- A + A1
      t_squared <- (nx * ny) / (nx + ny) * t(delta) %*% A %*% (delta)
      statistic <- t_squared * (nx + ny - i - 1) / (i * (nx + ny - 2))
      t_val[i] <- statistic
    }

    p_val_star[n, z] <- max(t_val)
  }
}


sums <- c()
for (i in 1:6) {
  sums[i] <- print(ApproxQuantile(hist(p_val_star[, i], breaks = 500, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95))
}

sums


###################### Brownian Error


p_val <- matrix(nrow = 10000, ncol = 6)
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)

for (z in 1:6) {
  for (n in 1:10000) {
    t <- seq(from = 0, to = 1, length.out = 25)
    x <- rbind(
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
    )

    y <- rbind(
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
    )

    nx <- nrow(x)
    ny <- nrow(y)
    delta <- colMeans(x) - colMeans(y)
    p <- ncol(x)
    Sx <- cov(x)
    Sy <- cov(y)
    S_pooled <- ((nx - 1) * Sx + (ny - 1) * Sy) / (nx + ny - 2)

    L <- eigen(S_pooled)$values
    V <- eigen(S_pooled)$vectors
    A <- matrix(rep(0), nrow = 25, ncol = 25)

    val <- c()
    for (i in 1:25) {
      val[i] <- L[i] / L[1]
    }

    k <- tail(which(val > 10^-16), 1)

    t_val <- c()
    for (i in 1:k) {
      L_prime <- as.numeric(solve(L[i]))
      A1 <- L_prime * V[, i] %*% t(V[, i])
      A <- A + A1
      t_squared <- (nx * ny) / (nx + ny) * t(delta) %*% A %*% (delta)
      statistic <- t_squared * (nx + ny - i - 1) / (i * (nx + ny - 2))
      t_val[i] <- statistic
    }

    p_val[n, z] <- max(t_val)
  }
}



sums <- c()
for (i in 1:6) {
  sums[i] <- print(ApproxQuantile(hist(p_val[, i], breaks = 500, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95))
}
