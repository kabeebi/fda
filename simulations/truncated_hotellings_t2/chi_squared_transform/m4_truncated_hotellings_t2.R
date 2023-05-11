##################### independent errors

sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45
t_val_star <- matrix(nrow = 1000, ncol = 6)
t <- seq(from = 0, to = 1, length.out = 25)


for (z in 1:6) {
  for (n in 1:1000) {
    x <- rbind(
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
    )

    y <- rbind(
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
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

    k <- tail(which(val > 10^-13), 1)

    t_val <- c()
    for (i in 1:k) {
      L_prime <- as.numeric(solve(L[i]))
      A1 <- L_prime * V[, i] %*% t(V[, i])
      A <- A + A1
      t_squared <- (nx * ny) / (nx + ny) * t(delta) %*% A %*% (delta)
      statistic <- (t_squared - i) * 1 / (sqrt(2 * i))
      t_val[i] <- statistic
    }
    t_val_star[n, z] <- max(t_val)
  }
}


sums <- c()
for (i in 1:6) {
  sums[i] <- sum(t_val_star[, i] > 17.5) / 1000
}
print(sums)



###################### Brownian Error


p_val <- matrix(nrow = 1000, ncol = 6)
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)
t <- seq(from = 0, to = 1, length.out = 25)

for (z in 1:6) {
  for (n in 1:1000) {
    x <- rbind(
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
    )


    y <- rbind(
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[z]
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

    k <- tail(which(val > 10^-13), 1)

    t_val <- c()
    for (i in 1:k) {
      L_prime <- as.numeric(solve(L[i]))
      A1 <- L_prime * V[, i] %*% t(V[, i])
      A <- A + A1
      t_squared <- (nx * ny) / (nx + ny) * t(delta) %*% A %*% (delta)
      statistic <- (t_squared - i) * 1 / (sqrt(2 * i))
      t_val[i] <- statistic
    }

    p_val[n, z] <- max(t_val)
  }
}


sums <- c()
for (i in 1:6) {
  sums[i] <- sum(p_val[, i] > 17.5) / 1000
}
print(sums)

################ Graphs

par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
plot(c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25),
  results_independent[, 2],
  xlim = c(0, 0.2), ylim = c(0, 1), lty = 1, type = "o", col = "black", xlab = "dispersion parameter", ylab = "p-value"
)
lines(c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25), results_independent[, 4], lty = 1, type = "o", col = "blue")
lines(c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25), results_independent[, 6], lty = 1, type = "o", col = "red")
lines(c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25), results_independent[, 8], lty = 1, type = "o", col = "green")
legend(
  x = "bottomright", # Position
  legend = c("M1", "M2", "M3", "M4"), # Legend texts
  lty = c(1), # Line types
  col = c("black", "blue", "red", "green"), # Line colors
  lwd = 2
)


# to transform to table in r
print(xtable(results_independent))
