############################ Null Distribution

################################### Independent Error

final <- c()
t <- seq(0, 1, length.out = 25) # time points
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45
sums <- c()

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

    statistic <- c()

    for (i in 1:25) {
      statistic[i] <- t.test(x[, i], y[, i])$statistic
      max_t <- max(statistic[2:25])
    }

    final[n1] <- max_t
  }

  sums[j] <- print(ApproxQuantile(hist(final, breaks = 500, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95))
}











############################ Null Distribution
final <- c()
t <- seq(0, 1, length.out = 25) # time points
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)
sums <- c()

for (j in 1:6) {
  for (n1 in 1:10000) {
    X <- rbind(
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    Y <- rbind(
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t * (1 - t) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    statistic <- c()

    for (i in 1:25) {
      statistic[i] <- t.test(X[, i], Y[, i])$statistic
      max_t <- max(statistic[2:25])
    }

    final[n1] <- max_t
  }

  sums[j] <- print(ApproxQuantile(hist(final, breaks = 500, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95))
}
