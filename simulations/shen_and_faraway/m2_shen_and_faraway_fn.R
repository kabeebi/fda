##################### independent errors

t_val_star <- matrix(nrow = 1000, ncol = 6)
t <- seq(from = 0, to = 1, length.out = 25)

sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45

for (j in 1:6) {
  for (n in 1:1000) {
    x <- cbind(
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

    y <- cbind(
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

    data <- cbind(x, y)
    data <- as.matrix(data)
    label <- rep(1:2, each = 10)

    fanova1 <- fanova.tests(data, group.label = label, test = "FN")
    t_val_star[n, j] <- as.numeric(as.character(unlist(fanova1$FN[2])))
  }
}

sums <- c()
for (i in 1:6) {
  sums[i] <- sum(t_val_star[, i] < 0.05) / 1000
}
sums


###################### Brownian Error

t_val_star <- matrix(nrow = 1000, ncol = 6)
p_val <- matrix(nrow = 1000, ncol = 6)
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)

for (j in 1:6) {
  for (n in 1:1000) {
    x <- cbind(
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

    y <- cbind(
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

    data <- cbind(x, y)
    data <- as.matrix(data)
    label <- rep(1:2, each = 10)

    fanova1 <- fanova.tests(data, group.label = label, test = "FN")
    t_val_star[n, j] <- as.numeric(as.character(unlist(fanova1$FN[2])))
  }
}

sums <- c()
for (i in 1:6) {
  sums[i] <- sum(t_val_star[, i] < 0.05) / 1000
}
sums
