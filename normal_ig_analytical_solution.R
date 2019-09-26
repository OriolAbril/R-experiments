source("integrate_2D.R")

# use analytical solution
log_p_y_given_gam <- function(y, Z, V0, g, a, b) {
  n <- length(y)
  m <- nrow(V0)
  V <- t(Z) %*% Z + solve(V0) / g
  const <- (-n * log(pi) - m * log(g) + a * log(b)) / 2
  const_gam <- lgamma((n + a) / 2) - lgamma(a / 2)
  dets <- -(
    determinant(V0, logarithm = TRUE)$modulus + determinant(V, logarithm = TRUE)$modulus
  ) / 2
  last_term <- -(n + a) / 2 * log(b +
    t(y) %*% y -
    t(y) %*% Z %*% solve(V) %*% t(Z) %*% y)
  const + const_gam + dets + last_term
}

n <- 10
m <- 1
y <- matrix(rnorm(n))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
V0 <- diag(m)
g <- 2
a <- 6
b <- 4

log_p_y_given_gam(y, Z, V0, g, a, b)

# use numerical integration (in 2D) to check analytical result
log_joint_p_1d <- function(theta, phi, y, Z, V0, g, a, b) {
  n <- length(y)
  m <- 1
  Z <- drop(Z)
  V0 <- drop(V0)[1]
  const1 <- -(n + m) / 2 * log(2 * pi) - m / 2 * log(g)
  const2 <- (a / 2) * log(b / 2) - lgamma(a / 2)
  v0_det <- -determinant(V0, logarithm = TRUE)$modulus / 2
  phi_term <- -((n + m + a) / 2 + 1) * log(phi)
  yzt <- y - Z %*% theta
  last_term <- -(t(yzt) %*% yzt + t(theta) %*% solve(V0) %*% theta / g + b) / (2 * phi)
  const1 + const2 + v0_det + phi_term + last_term
}

integrate2(log_joint_p, lower = c(-Inf, 0), upper = c(Inf, Inf), y = y, Z = Z, V0 = V0, g = g, a = a, b = b)
