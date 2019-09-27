# initialization
source("normal_ig_functions.R")
library(mombf)
set.seed(1234)

# define data
n <- 100
x <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
x <- cbind(matrix(1, nrow = n, ncol = 1), x) # add constant term
m <- ncol(x)
theta <- matrix(c(0, 1, 1, 0), ncol = 1)
y <- x %*% theta + rnorm(n)
g <- 0.348
a <- .1
b <- .1

###########################
# Use analytical solution #
###########################
gammas <- gtools::permutations(n = 2, r = 4, v = c(FALSE, TRUE), repeats.allowed = TRUE)
covariates_count <- apply(gammas, 1, sum)
df <- data.frame(modelid = character(), pp_analytical = double(), stringsAsFactors = FALSE)
for (k in 2:m^2) {
  gamma_k <- gammas[k,]
  df[k, "modelid"] <- paste(c(1,2,3,4)[gamma_k], collapse=",")
  Z <- x[, gamma_k]
  V0 <- t(Z) %*% Z
  model_like <- log_p_y_given_gam(y, Z=Z, V0=V0, g=g, a=a, b=b)
  model_prior <- (extraDistr::dbbinom(sum(gamma_k), m^2, 1, 1, log=TRUE) -
                  log(sum(covariates_count == sum(gamma_k))))
  df[k, "pp_analytical"] <- model_like + model_prior
}
df <- df[2:m^2,]
df[, "pp_analytical"] <- exp(df[, "pp_analytical"] - max(df[, "pp_analytical"]))

###########################
#   Use mombf solution    #
###########################
# define model
# Default MOM prior on parameters
priorCoef <- zellnerprior(tau = g)
priorVar <- igprior(a,b)
# Beta-Binomial prior for model space
priorDelta <- modelbbprior(1, 1)


# Model selection
fit1 <- modelSelection(y = y, x = x, priorCoef = priorCoef, priorDelta = priorDelta,
                       priorVar = priorVar, center = FALSE, scale = FALSE)

# Posterior model probabilities
post <- postProb(fit1)
post <- post[order(as.numeric(rownames(post))), ]
df[, "pp_mombf"] <- exp(log(post[, "pp"]) - max(log(post[, "pp"])))[2:m^2]
df
