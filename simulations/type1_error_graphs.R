df1 <- read_xlsx("Results for the Simulation.xlsx")[5:39, 2:9]
df2 <- read_xlsx("Results for the Simulation.xlsx")[5:39, 11:18]


# Graph of independent errors

df1_M1 <- df1[c(2, 7, 12, 17, 22, 27, 32), ]
colnames(df1_M1) <- c("Test", "Model", 0.0225, 0.0450, 0.0675, 0.090, 0.180, 0.360)
df1_M1 <- as.data.frame(df1_M1)

df1_m <- melt(df1_M1, id.vars = c("Test", "Model"))
df1_m$variable <- as.numeric(levels(df1_m$variable))[df1_m$variable]

ggplot(df1_m, aes(x = variable, y = value, colour = Test)) +
  geom_point() +
  geom_line() +
  xlab("Standard Deviation") +
  ylab("Type 1 Error")


# Graph of brownian errors

df2_M1 <- df2[c(2, 7, 12, 17, 22, 27, 32), ]
colnames(df2_M1) <- c("Test", "Model", 0.05, 0.1, 0.15, 0.2, 0.4, 0.8)
df2_M1 <- as.data.frame(df2_M1)

df2_m <- melt(df2_M1, id.vars = c("Test", "Model"))
df2_m$variable <- as.numeric(levels(df2_m$variable))[df2_m$variable]

ggplot(df2_m, aes(x = variable, y = value, colour = Test)) +
  geom_point() +
  geom_line() +
  xlab("Standard Deviation") +
  ylab("Type 1 Error")
