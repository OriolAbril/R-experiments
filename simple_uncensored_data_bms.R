# Attempt to implement a simple Gibbs sampler to use in
# Bayesian variable selection in linear regression cases

# Define model parameters
# True parameters
n <- 1000
# r <- 3
gamma_true <- c(1, 0, 0, 1)
beta_true <- c(2, -5)
sigma_true <- 2
alpha_true <- beta_true / sigma_true
rho_true <- -log(sigma_true)

# Prior parameters
g_L <- .1
# g_S <- n / r
a_sigma <- 4
b_sigma <- 1

# Data generation
x <- 1:n
x0 <- rep(1, n)
x1 <- -3 * x - 2
x2 <- x^2
x3 <- -x^3 + 20
x4 <- log(x) + 1

X_orig <- cbind(x1, x2, x3, x4)
X_orig <- matrix(rnorm(n * 4), n, 4)
X <- t((X_orig - apply(X_orig, 1, mean)) / apply(X_orig, 1, sd))
X_orig <- t(X_orig)
X_gam <- X[as.logical(gamma_true), ]
y <- t(X_gam) %*% beta_true + rnorm(n, mean = 0, sd = sigma_true)

###### Function definitions #######
reparametrize <- function(x) {
  n_x <- length(x)
  sigma <- exp(-x[n_x])
  c(x[1:n_x - 1] / sigma, sigma)
}

# priors
ig_const <- a_sigma / 2 * log(b_sigma / 2) - log(gamma(a_sigma / 2)) + log(2)
log_prior_ig <- function(rho, a = a_sigma, b = b_sigma, const = ig_const) {
  var_term <- a * rho - .5 * b * exp(2 * rho)
  const + var_term
}

log_prior_normal <- function(x, X, scale_factor = g_L) {
  # x and X input must have already been masked, i.e. X=X[mask, ]
  X2 <- apply(X * X, 1, sum)
  const <- -length(x) / 2 * log(2 * pi * scale_factor)
  const_sum <- .5 * sum(log(X2))
  var_term <- -1 / (2 * scale_factor) * sum(x^2 * X2)
  const + const_sum + var_term
}

log_prior_zellner <- function(alpha, rho, X) {
  log_prior_ig(rho) + log_prior_normal(alpha, X)
}

# prior gradients
grad_log_prior_ig <- function(rho, a = a_sigma, b = b_sigma) {
  a - b * exp(2 * rho)
}

grad_log_prior_normal <- function(x, X, scale_factor = g_L) {
  X2 <- apply(X * X, 1, sum)
  -1 / scale_factor * x * X2
}

grad_log_prior_zellner <- function(alpha, rho, X) {
  c(grad_log_prior_normal(alpha, X), grad_log_prior_ig(rho))
}

# likelihood
like_const <- -n / 2 * log(2 * pi)
log_likelihood <- function(y, alpha, rho, X, const = like_const, n_ = n) {
  const + n_ * rho - .5 * sum((exp(rho) * y - t(X) %*% alpha)^2)
}

# likelihood gradients
grad_log_likelihood <- function(y, alpha, rho, X, n_ = n) {
  grad_alpha <- X %*% (exp(rho) * y - t(X) %*% alpha)
  grad_rho <- n_ - exp(rho) * sum(y * (exp(rho) * y - t(X) %*% alpha))
  c(grad_alpha, grad_rho)
}

# combining everything together
log_posterior_fun <- function(theta, y, X) {
  alpha <- theta[1:length(theta) - 1]
  rho <- theta[length(theta)]
  result <- log_prior_zellner(alpha, rho, X) + log_likelihood(y, alpha, rho, X)
  -result
}

grad_log_posterior <- function(theta, y, X) {
  alpha <- theta[1:length(theta) - 1]
  rho <- theta[length(theta)]
  result <- grad_log_prior_zellner(alpha, rho, X)
  result <- result + grad_log_likelihood(y, alpha, rho, X)
  -result
}

true_pars <- c(alpha_true, rho_true)
true_pars

"like"
log_likelihood(y, c(0, 0), -3, X_gam)
"prior"
log_prior_zellner(c(0, 0), -3, X_gam)
log_prior_ig(-3)
log_prior_normal(c(0, 0), X_gam)
"posts"
log_posterior_fun(true_pars, y, X_gam)
log_posterior_fun(c(0, 0, -3), y, X_gam)
log_posterior_fun(true_pars + 1, y, X_gam)

reparametrize(
  optim(true_pars, log_posterior_fun, grad_log_posterior, y, X_gam)$par
)
reparametrize(
  optim(c(2, -5, 1), log_posterior_fun,
    grad_log_posterior, y, X_gam,
    method = "BFGS"
  )$par
)
reparametrize(
  optim(c(2, 1, -2, -5, 1), log_posterior_fun, grad_log_posterior, y, X)$par
)
reparametrize(
  optim(c(2, 1, -2, -5, 1), log_posterior_fun,
    grad_log_posterior, y, X,
    method = "BFGS"
  )$par
)

# Gibbs sampling
