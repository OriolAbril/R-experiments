paraboloid <- function(x) {
  2 * (x[1] - 3) ^ 2 + (x[2] - 5) ^ 2
}

paraboloid_grad <- function(x) {
  c(4 * (x[1] - 3), 2 * (x[2] - 5))
}

optim(c(0, 0), paraboloid)
optim(c(100, -345), paraboloid, paraboloid_grad, method="BFGS")
