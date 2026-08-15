## Spatial dependence
x1 = 1:10
x2 = 1:10

don1 <- data2(n1 = 10, n2 = 10, a = 0.1, sigma = 0.1)
don2 <- data2(n1 = 10, n2 = 10, a = 200, sigma = 0.1)

# Extract the column to visualize
mat1 <- matrix(don1[, 1], 10, 10)
mat2 <- matrix(don2[, 1], 10, 10)

# Display
windows(width = 10, height = 10)
par(mfrow = c(1, 2))  # Two plots side by side

image.plot(x1, x2, mat1, main = "a = 0.1", zlim = c(-3, 4))
image.plot(x1, x2, mat2, main = "a = 200", zlim = c(-3, 4))


# Smoothing of functional data
windows()

a = seq(0, 1, 1/364)

don3 <- data2(
  n1 = 10,
  n2 = 10,
  a = 0.1,
  sigma = 0.1
)[, 1:365]

matplot(
  a,
  t(don3),
  type = "l",
  lty = 1,
  col = "black",
  xlab = "Time",
  ylab = "X(t)",
  main = "Functional trajectory of X(t)"
)

# mat5 = as.matrix(don1)
# b <- array(don1, c(n1 = 10, n2 = 10))

plot(
  don1,
  type = "l"
)


## Response variable: replace your Z with the substituted variable
don <- as.matrix(data2(n1, n2, a, sigma))

x1 = 1:10
x2 = 1:10

don4 <- data2(
  n1 = 10,
  n2 = 10,
  a = 0.1,
  sigma = 0.1
)[, p + 3]

don5 <- data2(
  n1 = 10,
  n2 = 10,
  a = 0.1,
  sigma = 5
)[, p + 3]

mat4 <- matrix(don4, 10, 10)
mat5 <- matrix(don5, 10, 10)

windows(width = 10, height = 10)
par(mfrow = c(1, 2))  # Two plots side by side

image.plot(
  x1,
  x2,
  mat4,
  main = "Response variable when a = 0.1 and sigma = 0.1",
  zlim = c(-9, 30)
)

image.plot(
  x1,
  x2,
  mat5,
  main = "Response variable when a = 0.1 and sigma = 5",
  zlim = c(-9, 30)
)
