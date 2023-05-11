sigma_star <- c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25)
z <- 2
t1 <- c()

for (j in 1:100000) {
  t <- seq(from = 0, to = 1, length.out = 25)
  x <- rbind(
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z])
  )
  colnames(x) <- t

  y <- rbind(
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z]),
    t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma_star[z])
  )
  colnames(y) <- t

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

  t1[j] <- max(t_val)
}

hist(t1, breaks = 100, col = "lightblue", main = "Null Distribution", xlab = "Truncated Hotellings T^2 with F-Transform", ylab = "Frequency", xlim = c(0, 3))
ApproxQuantile(hist(t1, breaks = 100, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95)
