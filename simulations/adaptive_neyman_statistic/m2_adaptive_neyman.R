

################################### M1
################################### Independent Error

sigma_star <- c(0.2/25, 1/25, 1.8/25, 2.6/25, 3.4/25, 4.2/25, 5/25)

t <- seq(from = 0, to = 1, length.out = 25)
n_1 <- 10
n_2 <- 10
n <- 10
t_an_star <- matrix(nrow = 400, ncol = 7)

for (z in 1:7) {

for (j in 1:400) {
   x <- rbind(t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^1*(1-t)^(6-1) + rnorm(25, mean=0, sd=sigma_star[z]) )
    colnames(x) <- t
    
    y <- rbind(t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]),
               t^2*(1-t)^(6-2) + rnorm(25, mean=0, sd=sigma_star[z]) )
    colnames(y) <- t

x_trans <- fft(x)
y_trans <- fft(y)

x <- cbind ( c(rbind( Re(x_trans[1,]), Im(x_trans[1,])))[-2],
             c(rbind( Re(x_trans[2,]), Im(x_trans[2,])))[-2],
             c(rbind( Re(x_trans[3,]), Im(x_trans[3,])))[-2],
             c(rbind( Re(x_trans[4,]), Im(x_trans[4,])))[-2],
             c(rbind( Re(x_trans[5,]), Im(x_trans[5,])))[-2],
             c(rbind( Re(x_trans[6,]), Im(x_trans[6,])))[-2],
             c(rbind( Re(x_trans[7,]), Im(x_trans[7,])))[-2],
             c(rbind( Re(x_trans[8,]), Im(x_trans[8,])))[-2],
             c(rbind( Re(x_trans[9,]), Im(x_trans[9,])))[-2], 
             c(rbind( Re(x_trans[10,]), Im(x_trans[10,])))[-2] )
x <- x[c(1:25),]
x <- t(x)

y <- cbind ( c(rbind( Re(y_trans[1,]), Im(y_trans[1,])))[-2],
             c(rbind( Re(y_trans[2,]), Im(y_trans[2,])))[-2],
             c(rbind( Re(y_trans[3,]), Im(y_trans[3,])))[-2],
             c(rbind( Re(y_trans[4,]), Im(y_trans[4,])))[-2],
             c(rbind( Re(y_trans[5,]), Im(y_trans[5,])))[-2],
             c(rbind( Re(y_trans[6,]), Im(y_trans[6,])))[-2],
             c(rbind( Re(y_trans[7,]), Im(y_trans[7,])))[-2],
             c(rbind( Re(y_trans[8,]), Im(y_trans[8,])))[-2],
             c(rbind( Re(y_trans[9,]), Im(y_trans[9,])))[-2], 
             c(rbind( Re(y_trans[10,]), Im(y_trans[10,])))[-2] )
y <- y[c(1:25),]
y <- t(y)

xbar <- c()
for (i in 1:25) {
  xbar[i] <- mean(x[,i])
}

ybar <- c()
for (i in 1:25) {
  ybar[i] <- mean(y[,i])
}

sigma1 <- c()
for (i in 1:25) {
  sigma1[i] <- var(x[,i])
}

sigma2 <- c()
for (i in 1:25) {
  sigma2[i] <- var(y[,i])
}

#standardized difference
z_sd <- ( ( sigma1/n_1 ) + (sigma2/n_2) )^(-1/2) * (xbar - ybar)

#standardized difference transformed
t_star <- 1/(sqrt(2*25)) * sum(z_sd^2 -1)

t_an_star[j,z] <- sqrt(2*log(log(n)))*t_star - ( (2*log(log(n))  + 0.5 *log(log(log(n))) - 0.5*log(4*pi)) )

}

}
colnames(t_an_star) <- sigma_star


################################### Brownian Error

sigma <- c(0.2, 1, 1.8, 2.6, 3.4, 4.2, 5)

t <- seq(from = 0, to = 1, length.out = 25)
n_1 <- 10
n_2 <- 10
n <- 10
t_an <- matrix(nrow = 400, ncol = 7)

for (z in 1:7) {

for (j in 1:400) {
    x <- rbind(t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])),
               t^1*(1-t)^(6-1) + cumsum(rnorm(25,0,sigma[z])) )
    colnames(x) <- t
    
    y <- rbind(t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])),
               t^2*(1-t)^(6-2) + cumsum(rnorm(25,0,sigma[z])) )
    colnames(y) <- t


x_trans <- fft(x)
y_trans <- fft(y)

x <- cbind ( c(rbind( Re(x_trans[1,]), Im(x_trans[1,])))[-2],
             c(rbind( Re(x_trans[2,]), Im(x_trans[2,])))[-2],
             c(rbind( Re(x_trans[3,]), Im(x_trans[3,])))[-2],
             c(rbind( Re(x_trans[4,]), Im(x_trans[4,])))[-2],
             c(rbind( Re(x_trans[5,]), Im(x_trans[5,])))[-2],
             c(rbind( Re(x_trans[6,]), Im(x_trans[6,])))[-2],
             c(rbind( Re(x_trans[7,]), Im(x_trans[7,])))[-2],
             c(rbind( Re(x_trans[8,]), Im(x_trans[8,])))[-2],
             c(rbind( Re(x_trans[9,]), Im(x_trans[9,])))[-2], 
             c(rbind( Re(x_trans[10,]), Im(x_trans[10,])))[-2] )
x <- x[c(1:25),]
x <- t(x)

y <- cbind ( c(rbind( Re(y_trans[1,]), Im(y_trans[1,])))[-2],
             c(rbind( Re(y_trans[2,]), Im(y_trans[2,])))[-2],
             c(rbind( Re(y_trans[3,]), Im(y_trans[3,])))[-2],
             c(rbind( Re(y_trans[4,]), Im(y_trans[4,])))[-2],
             c(rbind( Re(y_trans[5,]), Im(y_trans[5,])))[-2],
             c(rbind( Re(y_trans[6,]), Im(y_trans[6,])))[-2],
             c(rbind( Re(y_trans[7,]), Im(y_trans[7,])))[-2],
             c(rbind( Re(y_trans[8,]), Im(y_trans[8,])))[-2],
             c(rbind( Re(y_trans[9,]), Im(y_trans[9,])))[-2], 
             c(rbind( Re(y_trans[10,]), Im(y_trans[10,])))[-2] )
y <- y[c(1:25),]
y <- t(y)

xbar <- c()
for (i in 1:25) {
  xbar[i] <- mean(x[,i])
}

ybar <- c()
for (i in 1:25) {
  ybar[i] <- mean(y[,i])
}

sigma1 <- c()
for (i in 1:25) {
  sigma1[i] <- var(x[,i])
}

sigma2 <- c()
for (i in 1:25) {
  sigma2[i] <- var(y[,i])
}

#standardized difference
z_sd <- ( ( sigma1/n_1 ) + (sigma2/n_2) )^(-1/2) * (xbar - ybar)

#standardized difference transformed
t_star <- 1/(sqrt(2*30)) * sum(z_sd^2 -1)

t_an[j,z] <- sqrt(2*log(log(n)))*t_star - ( (2*log(log(n))  + 0.5 *log(log(log(n))) - 0.5*log(4*pi)) )

}

}
colnames(t_an_star) <- sigma_star

colnames(t_an) <- sigma

#########################################################M1


M2 <- cbind(t_an_star,t_an)

#########################################################Finding the p-value and the acceptance rate in all the models

#Independent Errors

sampled <- cbind(sample(M2[, 1], size = 100, replace = FALSE),
                 sample(M2[, 2], size = 100, replace = FALSE),
                 sample(M2[, 3], size = 100, replace = FALSE),
                 sample(M2[, 4], size = 100, replace = FALSE),
                 sample(M2[, 5], size = 100, replace = FALSE),
                 sample(M2[, 6], size = 100, replace = FALSE),
                 sample(M2[, 7], size = 100, replace = FALSE) )

for (i in 1:7) {
  results_independent[3,i] <- mean(sampled[,i])
  results_independent[4,i] <- sum(sampled[,i] < 3.77)/100
}


#Brownian Errors

sampled <- cbind(sample(M2[, 8], size = 100, replace = FALSE),
                 sample(M2[, 9], size = 100, replace = FALSE),
                 sample(M2[, 10], size = 100, replace = FALSE),
                 sample(M2[, 11], size = 100, replace = FALSE),
                 sample(M2[, 12], size = 100, replace = FALSE),
                 sample(M2[, 13], size = 100, replace = FALSE),
                 sample(M2[, 14], size = 100, replace = FALSE))

for (i in 1:7) {
  results_brownian[3,i] <- mean(sampled[,i])
  results_brownian[4,i] <- sum(sampled[,i] < 3.77)/100
}



