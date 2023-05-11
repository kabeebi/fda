##################### Brownian Errors

final1 <- c()
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8)
sums <- c()
t <- seq(0, 1, length.out = 25) # time points

for (j in 1:6) {
  for (n2 in 1:1000) {
    # generate Brownian motion
    X <- rbind(
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      1 + (1 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    Y <- rbind(
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j],
      2 + (2 / 50) + BM(x = 0, t0 = 0, T = 1, N = 24) * sigma[j]
    )

    statistic <- c()

    for (i in 2:25) {
      statistic[i] <- t.test(X[, i], Y[, i])$statistic
      max_t <- max(statistic[2:25])
    }

    final1[n2] <- max_t
  }

  sums[j] <- sum(final1 > 2.73) / 1000
}


sums





################# independent errors

maximum_t_star <- c()
t <- seq(from = 0, to = 1, length.out = 25)
sigma <- c(0.05, 0.1, 0.15, 0.2, 0.4, 0.8) * 0.45
sums <- c()

for (j in 1:6) {
  for (n in 1:1000) {
    x <- rbind(
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      1 + (1 / 50) + rnorm(25, mean = 0, sd = sigma[j])
    )

    y <- rbind(
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j]),
      2 + (2 / 50) + rnorm(25, mean = 0, sd = sigma[j])
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
