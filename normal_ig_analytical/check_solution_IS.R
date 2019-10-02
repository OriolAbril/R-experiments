library("rstan")
source("../integrate_2D.R")
source("normal_ig_functions.R")

set.seed(14)

# define parameters used -----------------
n <- 100
m <- 5
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
V0 <- diag(m)
g <- .1
a <- 1
b <- 1

# use analytical solution ----------------
exp(log_p_y_given_gam(y, Z, V0, g, a, b)[1])

# use numerical approach -----------------
# use importance sampling
# get posterior samples of phi and theta with stan
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
stan_dat <- list(
  n = n, m = m, y = y[, 1], Z = Z, V0 = V0, g = g, a = a, b = b
)

fit <- purrr::quietly(rstan::stan)(
  file = "linreg_ig_normal.stan", data = stan_dat, verbose = FALSE
)$result
posterior <- rstan::extract(fit)

# use posterior samples to calculate integral
# with importance sampling
log_jp <- vector() # store log_joint_probability
log_pp <- vector() # store log_posterior of theta and phi
for (i in 1:dim(posterior$phi)) {
  theta_i <- posterior$theta[i, ]
  phi_i <- posterior$phi[i]
  log_jp[i] <- log_joint_p_1d(theta_i, phi_i, y, Z, V0, g, a, b)
  log_pp[i] <- log_post_theta_phi(theta_i, phi_i, y, Z, V0, g, a, b)
}
is_weights <- log_jp - log_pp - log(length(log_jp))

exp(matrixStats::logSumExp(is_weights))
