################################### M1
################################### Independent Error

sigma_star <- c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25)

t <- seq(from = 0, to = 1, length.out = 25)
n_1 <- 10
n_2 <- 10
n <- 20
t_an_star <- matrix(nrow = 400, ncol = 7)

for (z in 1:7) {
  for (j in 1:400) {
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

    # standardized difference transformed
    t_star <- 1 / (sqrt(2 * length(z_star))) * sum(z_star^2 - 1)

    t_an_star[j, z] <- max(t_star)
  }
}
colnames(t_an_star) <- sigma_star


################################### Brownian Error

sigma <- c(0.2, 1, 1.8, 2.6, 3.4, 4.2, 5)

t <- seq(from = 0, to = 1, length.out = 25)
n_1 <- 10
n_2 <- 10
n <- 25
t_an <- matrix(nrow = 400, ncol = 7)

for (z in 1:7) {
  for (j in 1:400) {
    x <- rbind(
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z]))
    )
    colnames(x) <- t

    y <- rbind(
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z])),
      t * (1 - t) + cumsum(rnorm(25, 0, sigma[z]))
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

    # standardized difference transformed
    t_star <- 1 / (sqrt(2 * length(z_star))) * sum(z_star^2 - 1)

    t_an_star[j, z] <- max(t_star)
  }
}
colnames(t_an_star) <- sigma_star

colnames(t_an) <- sigma


######################################################### M1


M1 <- cbind(t_an_star, t_an)

######################################################### Finding the p-value and the acceptance rate in all the models


sampled <- cbind(
  sample(M1[, 1], size = 100, replace = FALSE),
  sample(M1[, 2], size = 100, replace = FALSE),
  sample(M1[, 3], size = 100, replace = FALSE),
  sample(M1[, 4], size = 100, replace = FALSE),
  sample(M1[, 5], size = 100, replace = FALSE),
  sample(M1[, 6], size = 100, replace = FALSE),
  sample(M1[, 7], size = 100, replace = FALSE)
)

results_independent <- matrix(nrow = 8, ncol = 7)
colnames(results_independent) <- c(0.2 / 25, 1 / 25, 1.8 / 25, 2.6 / 25, 3.4 / 25, 4.2 / 25, 5 / 25)

for (i in 1:7) {
  results_independent[1, i] <- mean(sampled[, i])
  results_independent[2, i] <- sum(sampled[, i] < 3.77) / 100
}


# Brownian Errors

sampled <- cbind(
  sample(M1[, 8], size = 100, replace = FALSE),
  sample(M1[, 9], size = 100, replace = FALSE),
  sample(M1[, 10], size = 100, replace = FALSE),
  sample(M1[, 11], size = 100, replace = FALSE),
  sample(M1[, 12], size = 100, replace = FALSE),
  sample(M1[, 13], size = 100, replace = FALSE),
  sample(M1[, 14], size = 100, replace = FALSE)
)

results_brownian <- matrix(nrow = 8, ncol = 7)
colnames(results_brownian) <- c(0.2, 1, 1.8, 2.6, 3.4, 4.2, 5)

for (i in 1:7) {
  results_brownian[1, i] <- mean(sampled[, i])
  results_brownian[2, i] <- sum(sampled[, i] < 3.77) / 100
}
