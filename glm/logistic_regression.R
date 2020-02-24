library("RegressionFactory")
source("helper_functions.R")

loglike.logistic <- function(beta, X, y, fgh) {
  regfac.expand.1par(beta, X, y, fbase1.binomial.logit, fgh, n=1)
}

logprior.logistic <- function(beta, mu.beta, sd.beta, fgh) {
  f <- sum(dnorm(beta, mu.beta, sd.beta, log=TRUE))
  if (fgh==0) return (f)
  g <- -(beta-mu.beta)/sd.beta^2
  if (fgh==1) return (list(f=f, g=g))
  h <- diag(-1/sd.beta^2, nrow=length(beta))
  return (list(f=f, g=g, h=h))
}

logpost.logistic <- function(beta, X, y, mu.beta, sd.beta, fgh) {
  ret.loglike <- loglike.logistic(beta, X, y, fgh)
  ret.logprior <- logprior.logistic(beta, mu.beta, sd.beta, fgh)
  regfac.merge(ret.loglike, ret.logprior, fgh=fgh)
}

n_rep <- 10
optim_iter <- 2
N <- 1000
K <- 5
X <- matrix(runif(N*K, min=-0.5, max=+0.5), ncol=K)
beta <- c(1.45, 0, 0, -1.13, .85)
y <- rbinom(N, size = 1, prob = 1/(1+exp(-X%*%beta)))
beta.glm <- glm(y~X-1, family="binomial")$coefficients
beta.glm

verbose = FALSE
gammas <- list_vars(K)
for (idx in 1:n_rep) {
  min_lapprox_str <- ""
  min_lapprox <- 1e19
  par <- rnorm(length(beta), sd=1)
  cat(paste("beta initial guess:", paste(par, collapse=", "), "\n"))
  for (row in 2:nrow(gammas)) {
    gamma <- gammas[row, ]
    if (verbose) {cat(paste(gamma_to_str(gamma), "\t"))}
    X_j <- X[, gamma, drop=FALSE]
    par_j <- par[gamma, drop=FALSE]
    lapprox <- laplace_logistic_post(par_j, logpost.logistic, optim_iter, X_j, y, 0, 1000)
    if (!is.na(lapprox) & lapprox < min_lapprox) {
      min_lapprox_str <- gamma_to_str(gamma)
      min_lapprox <- lapprox
    }
    if (verbose) {cat(paste(lapprox), "\n")}
  }
  cat(paste(min_lapprox_str, "\n"))
}

