sigma_star <- c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25)

t <- seq(from = 0, to = 1, length.out = 25)
n_1 <- 10
n_2 <- 10
n <- 20
z <- 4
t2 <- c()
t3 <- c()

for (z in 1:7) {
  for (j in 1:100000) {
    x <- rbind(
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z])
    )
    colnames(x) <- t

    y <- rbind(
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z]),
      t * (1 - t) + rnorm(25, mean = 0, sd = sigma_star[z])
    )
    colnames(y) <- t

    x_trans <- fft(x)
    y_trans <- fft(y)

    x_trans_1 <- matrix(0, nrow = 19, ncol = 25)
    for (i in 1:25) {
      x_trans_1[, i] <- cbind(c(rbind(Re(x_trans[, i]), Im(x_trans[, i])))[-2])
    }
    x_trans_1 <- x_trans_1[c(1:10), ]

    y_trans_1 <- matrix(0, nrow = 19, ncol = 25)
    for (i in 1:25) {
      y_trans_1[, i] <- cbind(c(rbind(Re(y_trans[, i]), Im(y_trans[, i])))[-2])
    }
    y_trans_1 <- y_trans_1[c(1:10), ]


    xbar <- c()
    for (i in 1:25) {
      xbar[i] <- mean(x_trans_1[, i])
    }

    ybar <- c()
    for (i in 1:25) {
      ybar[i] <- mean(y_trans_1[, i])
    }

    sigma1 <- c()
    for (i in 1:25) {
      sigma1[i] <- var(x_trans_1[, i])
    }

    sigma2 <- c()
    for (i in 1:25) {
      sigma2[i] <- var(y_trans_1[, i])
    }

    # standardized difference
    z_sd <- ((n_1^(-1) * sigma1) + (n_2^(-1) * sigma2))^(-1 / 2) * (xbar - ybar)

    # find the value for k hat, z_star is the different values for test statistic under different k values
    z_star <- c()
    for (i in 1:25) {
      z_star[i] <- i^(-0.5) * sum(z_sd[i]^2 - 1)
    }

    z_star <- z_star[seq_len(which.max(z_star))]

    # standardized difference transformed
    t_star <- 1 / (sqrt(2 * length(z_star))) * sum(z_star^2 - 1)

    t2[j] <- max(t_star)
  }

  t3[z] <- ApproxQuantile(hist(t2, breaks = 100, col = "lightblue", main = "Null Distribution", xlab = "Maximum t-statistic", ylab = "Frequency"), 0.95)
}


# after running the whole thing for 5 different
# 0.2/25: -0.25
# 1/25: 1.5
# 1.8/25: 1.5
# 2.6/25,
# 3.4/25,
# 4.2/25,
# 5/25
