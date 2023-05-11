############################ 3
# Single Adaptive Neyman

x1 <- rnorm(15, mean = 3, sd = 1)
x2 <- rnorm(15, mean = 3, sd = 1)
x3 <- rnorm(15, mean = 3, sd = 1)
x4 <- rnorm(15, mean = 3, sd = 1)
x5 <- rnorm(15, mean = 3, sd = 1)

y1 <- rnorm(15, mean = 3, sd = 1)
y2 <- rnorm(15, mean = 3, sd = 1)
y3 <- rnorm(15, mean = 3, sd = 1)
y4 <- rnorm(15, mean = 3, sd = 1)
y5 <- rnorm(15, mean = 3, sd = 1)

x_trans <- cbind(fft(x1), fft(x2), fft(x3), fft(x4), fft(x5))
y_trans <- cbind(fft(y1), fft(y2), fft(y3), fft(y4), fft(y5))

x <- cbind(
  c(rbind(Re(x_trans[, 1]), Im(x_trans[, 1])))[-2],
  c(rbind(Re(x_trans[, 2]), Im(x_trans[, 2])))[-2],
  c(rbind(Re(x_trans[, 3]), Im(x_trans[, 3])))[-2],
  c(rbind(Re(x_trans[, 4]), Im(x_trans[, 4])))[-2],
  c(rbind(Re(x_trans[, 5]), Im(x_trans[, 5])))[-2]
)
x <- x[c(1:15), ]

y <- cbind(
  c(rbind(Re(y_trans[, 1]), Im(y_trans[, 1])))[-2],
  c(rbind(Re(y_trans[, 2]), Im(y_trans[, 2])))[-2],
  c(rbind(Re(y_trans[, 3]), Im(y_trans[, 3])))[-2],
  c(rbind(Re(y_trans[, 4]), Im(y_trans[, 4])))[-2],
  c(rbind(Re(y_trans[, 5]), Im(y_trans[, 5])))[-2]
)
y <- y[c(1:15), ]

xbar <- c()
for (i in 1:15) {
  xbar[i] <- mean(x[i, ])
}

ybar <- c()
for (i in 1:15) {
  ybar[i] <- mean(y[i, ])
}

sigma1 <- c()
for (i in 1:15) {
  sigma1[i] <- var(x[i, ])
}

sigma2 <- c()
for (i in 1:15) {
  sigma2[i] <- var(y[i, ])
}

z <- (sqrt((1 / 5) * sigma1 + (1 / 5) * sigma2))^(-1) * (xbar - ybar)

t <- 1 / (sqrt(2 * 15)) * sum(z^2 - 1)





#############################
# Simulation

t_values <- c()

for (n in 1:10000) {
  x1 <- rnorm(15, mean = 3, sd = 1)
  x2 <- rnorm(15, mean = 3, sd = 1)
  x3 <- rnorm(15, mean = 3, sd = 1)
  x4 <- rnorm(15, mean = 3, sd = 1)
  x5 <- rnorm(15, mean = 3, sd = 1)

  y1 <- rnorm(15, mean = 3, sd = 1)
  y2 <- rnorm(15, mean = 3, sd = 1)
  y3 <- rnorm(15, mean = 3, sd = 1)
  y4 <- rnorm(15, mean = 3, sd = 1)
  y5 <- rnorm(15, mean = 3, sd = 1)

  x_trans <- cbind(fft(x1), fft(x2), fft(x3), fft(x4), fft(x5))
  y_trans <- cbind(fft(y1), fft(y2), fft(y3), fft(y4), fft(y5))

  x <- cbind(
    c(rbind(Re(x_trans[, 1]), Im(x_trans[, 1])))[-2],
    c(rbind(Re(x_trans[, 2]), Im(x_trans[, 2])))[-2],
    c(rbind(Re(x_trans[, 3]), Im(x_trans[, 3])))[-2],
    c(rbind(Re(x_trans[, 4]), Im(x_trans[, 4])))[-2],
    c(rbind(Re(x_trans[, 5]), Im(x_trans[, 5])))[-2]
  )
  x <- x[c(1:15), ]

  y <- cbind(
    c(rbind(Re(y_trans[, 1]), Im(y_trans[, 1])))[-2],
    c(rbind(Re(y_trans[, 2]), Im(y_trans[, 2])))[-2],
    c(rbind(Re(y_trans[, 3]), Im(y_trans[, 3])))[-2],
    c(rbind(Re(y_trans[, 4]), Im(y_trans[, 4])))[-2],
    c(rbind(Re(y_trans[, 5]), Im(y_trans[, 5])))[-2]
  )
  y <- y[c(1:15), ]

  xbar <- c()
  for (i in 1:15) {
    xbar[i] <- mean(x[i, ])
  }

  ybar <- c()
  for (i in 1:15) {
    ybar[i] <- mean(y[i, ])
  }

  sigma1 <- c()
  for (i in 1:15) {
    sigma1[i] <- var(x[i, ])
  }

  sigma2 <- c()
  for (i in 1:15) {
    sigma2[i] <- var(y[i, ])
  }

  z <- (sqrt((1 / 5) * sigma1 + (1 / 5) * sigma2))^(-1) * (xbar - ybar)
  t <- 1 / (sqrt(2 * 15)) * sum(z^2 - 1)
  t_values[n] <- t
}

count(t_values > 2.97)






###############################

t_values <- c()

for (n in 1:10000) {
  x1 <- rnorm(15, mean = 3, sd = 1)
  x2 <- rnorm(15, mean = 3, sd = 1)
  x3 <- rnorm(15, mean = 3, sd = 1)
  x4 <- rnorm(15, mean = 3, sd = 1)
  x5 <- rnorm(15, mean = 3, sd = 1)

  y1 <- rnorm(15, mean = 3, sd = 1)
  y2 <- rnorm(15, mean = 3, sd = 1)
  y3 <- rnorm(15, mean = 3, sd = 1)
  y4 <- rnorm(15, mean = 3, sd = 1)
  y5 <- rnorm(15, mean = 3, sd = 1)

  x_trans <- cbind(fft(x1), fft(x2), fft(x3), fft(x4), fft(x5))
  y_trans <- cbind(fft(y1), fft(y2), fft(y3), fft(y4), fft(y5))

  x <- cbind(
    c(rbind(Re(x_trans[, 1]), Im(x_trans[, 1])))[-2],
    c(rbind(Re(x_trans[, 2]), Im(x_trans[, 2])))[-2],
    c(rbind(Re(x_trans[, 3]), Im(x_trans[, 3])))[-2],
    c(rbind(Re(x_trans[, 4]), Im(x_trans[, 4])))[-2],
    c(rbind(Re(x_trans[, 5]), Im(x_trans[, 5])))[-2]
  )
  x <- x[c(1:15), ]

  y <- cbind(
    c(rbind(Re(y_trans[, 1]), Im(y_trans[, 1])))[-2],
    c(rbind(Re(y_trans[, 2]), Im(y_trans[, 2])))[-2],
    c(rbind(Re(y_trans[, 3]), Im(y_trans[, 3])))[-2],
    c(rbind(Re(y_trans[, 4]), Im(y_trans[, 4])))[-2],
    c(rbind(Re(y_trans[, 5]), Im(y_trans[, 5])))[-2]
  )
  y <- y[c(1:15), ]

  xbar <- c()
  for (i in 1:15) {
    xbar[i] <- mean(x[i, ])
  }

  ybar <- c()
  for (i in 1:15) {
    ybar[i] <- mean(y[i, ])
  }

  sigma1 <- c()
  for (i in 1:15) {
    sigma1[i] <- var(x[i, ])
  }

  sigma2 <- c()
  for (i in 1:15) {
    sigma2[i] <- var(y[i, ])
  }

  z <- (sqrt((1 / 5) * sigma1 + (1 / 5) * sigma2))^(-1) * (xbar - ybar)
  t <- 1 / (sqrt(2 * 15)) * sum(z^2 - 1)
  t_values[n] <- t
}

count(t_values > 2.97)

hist(t_values)



##################################
# Corrected simulation for Adaptive Neyman:

t_values <- c()

x1 <- rnorm(15, mean = 3, sd = 1)
x2 <- rnorm(15, mean = 3, sd = 1)
x3 <- rnorm(15, mean = 3, sd = 1)
x4 <- rnorm(15, mean = 3, sd = 1)
x5 <- rnorm(15, mean = 3, sd = 1)

y1 <- rnorm(15, mean = 3, sd = 1)
y2 <- rnorm(15, mean = 3, sd = 1)
y3 <- rnorm(15, mean = 3, sd = 1)
y4 <- rnorm(15, mean = 3, sd = 1)
y5 <- rnorm(15, mean = 3, sd = 1)

total <- cbind(x1, x2, x3, x4, x5, y1, y2, y3, y4, y5)

for (n in 1:10000) {
  x <- total[, sample(ncol(total), size = 5), drop = FALSE]
  y <- total[, sample(ncol(total), size = 5), drop = FALSE]

  x_trans <- cbind(fft(x[, 1]), fft(x[, 2]), fft(x[, 3]), fft(x[, 4]), fft(x[, 5]))
  y_trans <- cbind(fft(y[, 1]), fft(y[, 2]), fft(y[, 3]), fft(y[, 4]), fft(y[, 5]))

  x <- cbind(
    c(rbind(Re(x_trans[, 1]), Im(x_trans[, 1])))[-2],
    c(rbind(Re(x_trans[, 2]), Im(x_trans[, 2])))[-2],
    c(rbind(Re(x_trans[, 3]), Im(x_trans[, 3])))[-2],
    c(rbind(Re(x_trans[, 4]), Im(x_trans[, 4])))[-2],
    c(rbind(Re(x_trans[, 5]), Im(x_trans[, 5])))[-2]
  )
  x <- x[c(1:15), ]

  y <- cbind(
    c(rbind(Re(y_trans[, 1]), Im(y_trans[, 1])))[-2],
    c(rbind(Re(y_trans[, 2]), Im(y_trans[, 2])))[-2],
    c(rbind(Re(y_trans[, 3]), Im(y_trans[, 3])))[-2],
    c(rbind(Re(y_trans[, 4]), Im(y_trans[, 4])))[-2],
    c(rbind(Re(y_trans[, 5]), Im(y_trans[, 5])))[-2]
  )
  y <- y[c(1:15), ]

  xbar <- c()
  for (i in 1:15) {
    xbar[i] <- mean(x[i, ])
  }

  ybar <- c()
  for (i in 1:15) {
    ybar[i] <- mean(y[i, ])
  }

  sigma1 <- c()
  for (i in 1:15) {
    sigma1[i] <- var(x[i, ])
  }

  sigma2 <- c()
  for (i in 1:15) {
    sigma2[i] <- var(y[i, ])
  }

  z <- (sqrt((1 / 5) * sigma1 + (1 / 5) * sigma2))^(-1) * (xbar - ybar)
  t <- 1 / (sqrt(2 * 15)) * sum(z^2 - 1)
  t_values[n] <- t
}

count(t_values > 2.97)
