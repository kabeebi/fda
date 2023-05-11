##################### independent errors

maximum_t_star <- c()
t <- seq(from = 0, to = 1, length.out = 25)
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45
sums <- c()

for (j in 1:6) {
  for (n in 1:1000) {
    x <- rbind(
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j]),
      t^1 * (1 - t)^(6 - 1) + rnorm(25, mean = 0, sd = sigma[j])
    )

    y <- rbind(
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j]),
      t^2 * (1 - t)^(6 - 2) + rnorm(25, mean = 0, sd = sigma[j])
    )

    statistic <- c()

    for (i in 1:25) {
      statistic[i] <- t.test(x[, i], y[, i])$statistic
      max_t <- max(statistic)
    }
    maximum_t_star[n] <- max_t
  }

  sums[j] <- sum(maximum_t_star > 3.27) / 1000
}



sums







###################### Brownian Error

final1 <- c()
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)
sums <- c()
t <- seq(0, 1, length.out = 25) # time points

for (j in 1:6) {
  for (n2 in 1:1000) {
    # generate Brownian motion
    X <- rbind(
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^1 * (1 - t)^(6 - 1) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    Y <- rbind(
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      t^2 * (1 - t)^(6 - 2) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    statistic <- c()

    for (i in 1:25) {
      statistic[i] <- t.test(X[, i], Y[, i])$statistic
      max_t <- max(statistic[2:25])
    }

    final1[n2] <- max_t
  }

  sums[j] <- sum(final1 > 2.71) / 1000
}

sums
