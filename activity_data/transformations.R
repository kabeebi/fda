data <- read_excel("500001_Summary by week.xlsx")

data$time <- c(1:169)
plot(data$time, data$`Step count`)

ggplot(data, aes(time, `Step count`)) +
  geom_point() +
  geom_line() +
  xlab("Time") +
  ylab("Number of Steps")


# fourier
data$`Step count` <- as.numeric(data$`Step count`)
y <- data[, 6]
x <- 1:169
colnames(y) <- c("Step_Count")
y <- as.numeric(unlist(y))
y_fft <- fft(y)

df <- data.frame(
  freq = seq_along(y_fft),
  amplitude = Mod(y_fft) # use Mod() to get the magnitude of each Fourier coefficient
)

ggplot(df, aes(x = freq, y = amplitude)) +
  geom_line() +
  labs(x = "Frequency", y = "Amplitude")


# splines
install.packages("splines")
library(splines)

wt <- modwt(x, n.levels = 2, boundary = "periodic")
data <- as.ts(data)
plot(wt)


knots <- 4
bs_basis <- bs(y, knots, intercept = T)
fit <- lm(y ~ bs_basis)
pred <- predict(fit, data[, 9])

plot(1:169, y)
lines(data[, 9], pred, col = "red")

plot(bs(y,
  df = NULL, knots = 3, degree = 4, intercept = FALSE,
  Boundary.knots = range(data[, 6])
))




x <- seq(1, 169)
y <- c(72, 6, 78, 22, 100, 0, 26, 38, 1348, 262, 1294, 1208, 6, 0, 608, 36, 28, 28, 1636, 772, 660, 240, 0, 206, 20, 0, 56, 0, 0, 0, 0, 0, 74, 1076, 1546, 0, 944, 158, 0, 0, 16, 16, 1078, 1176, 530, 396, 274, 148, 334, 312, 48, 30, 0, 92, 0, 0, 102, 1392, 1222, 516, 14, 442, 18, 210, 50, 24, 862, 1240, 74, 1184, 580, 396, 298, 166, 120, 86, 42, 24, 0, 50, 0, 0, 0, 120, 0, 64, 192, 124, 0, 0, 296, 22, 84, 0, 0, 102, 0, 0, 0, 0, 0, 0, 0, 0, 128, 370, 588, 272, 954, 510, 120, 46, 148, 74, 80, 294, 456, 108, 0, 96, 0, 36, 0, 0, 0, 0, 32, 0, 98, 802, 1368, 0, 2, 540, 46, 0, 0, 0, 1698, 1516, 1068, 672, 240, 178, 78, 34, 34, 28, 0, 0, 0, 0, 30, 120, 574, 176, 272, 190, 0, 44, 204, 34, 388, 104, 340, 204, 0, 28, 0) # nolint

# Use bs() to fit a B-spline
b_spline <- bs(x, y, degree = 3, knots = c(50, 100, 150), Boundary.knots = c(1, length(x)))

# Plot the original data and the B-spline
plot(x, y, pch = 16)
lines(x, predict(b_spline[, 1]), col = "red")
